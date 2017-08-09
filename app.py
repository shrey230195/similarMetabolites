import os
from flask import Flask, render_template, request, redirect, url_for, send_from_directory,jsonify
from werkzeug import secure_filename
import pandas as pd
import numpy as np
from scipy import stats
# Initialize the Flask application
app = Flask(__name__)

global uploadedFile
global df
global PC
global LPC
global Plasmalogen  
# This is the path to the upload directory
app.config['UPLOAD_FOLDER'] = 'uploads/'
# These are the extension that we are accepting to be uploaded
app.config['ALLOWED_EXTENSIONS'] = set(['csv', 'xlsx'])

# For a given file, return whether it's an allowed type or not
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']

def get_similar_metabolites(dataset):
    #create dataframe in such a way that new dataframe is
    #left with sample columns and Accepted_Compound_ID columns
    #which is basically column on metabolites   
    testDataset=dataset.drop(Plasmalogen.columns[[0,1]], axis=1)
    # #turn the tables i.e
    # #Now each entry in col[Accepted_Compound_ID] is now a separate column
    # #i.e each metabolite is a separate column now having values as information
    # #of 1050 samples
    transDataset = testDataset.transpose()
    #data cleaning : getting rid of index set as column after transpose
    #setting the first row as column headers
    transDataset.columns = transDataset.iloc[0]
    transDataset = transDataset.ix[1:]
    #this can be easily optimized with look_ahead or a backtrack dictionary
    similar_metabolites=[]
    for column in transDataset.columns:    
        for head in transDataset.columns:        
            if column!=head:
                z_stat, p_val = stats.ranksums(transDataset[column],transDataset[head])
                #if score is greater than 0.05 
                # => metabolites are statistically simillar 
                if p_val>0.05:
                    arr={
                        "metabolite1":column,
                        "metabolite2":head,
                        "score":p_val
                    }
                    if len(similar_metabolites):
                        for item in similar_metabolites:                        
                            if (item['metabolite1']!=column and item['metabolite2']!=head)and(item['metabolite2']!=column and item['metabolite1']!=head):                                      
                                similar_metabolites.append(arr)
                    else:
                        similar_metabolites.append(arr)
    return similar_metabolites

@app.route('/')
def index():
    return render_template('index.html')


# Route that will process the file upload
@app.route('/upload', methods=['POST'])
def upload():  
    global uploadedFile  
    # Get the name of the uploaded file
    file = request.files['file']
    # Check if the file is one of the allowed types/extensions
    if file and allowed_file(file.filename):
        # Make the filename safe, remove unsupported chars
        filename = secure_filename(file.filename)
        uploadedFile=filename
        # Move the file form the temporal folder to
        # the upload folder we setup
        file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
        # Redirect the user to the uploaded_file route, which        
        return redirect(url_for('uploaded_file'))

@app.route('/output')
def uploaded_file():
    global uploadedFile
    if(not uploadedFile):
       return redirect(url_for('index')) 
    return render_template('result.html')

@app.route('/processfile')
def process_file():
    global uploadedFile
    global df 
    global PC
    global LPC
    global Plasmalogen 
    if(not uploadedFile):
       return redirect(url_for('index')) 
    df = pd.read_excel('./uploads/'+uploadedFile, sep='none')
    #replace spaces with _ in column heads
    df.columns = [x.strip().replace(' ', '_') for x in df.columns]
    # data cleaning - get rid of noice
    df = df[df['Accepted_Compound_ID'].str.contains('_')==True]
    # ending with PC dataset
    PC = df[df['Accepted_Compound_ID'].str.endswith(' PC')==True]
    #ending with LPC dataset
    LPC = df[df['Accepted_Compound_ID'].str.endswith('LPC')==True]
    #ending with Plasmalogen dataset
    Plasmalogen = df[df['Accepted_Compound_ID'].str.endswith('plasmalogen')==True]
    #round off 'Retention_time_(min) to nearesr natural number
    value1 = df['Retention_time_(min)'].round()
    #insert Retention_time_round_off(min) at 3rd postion adjescent to Retention_time_(min)
    df.insert(2, 'Retention_time_round_off(min)', value1)    

    #drop categorical data for more accurate mean
    PC1=PC.drop(PC.columns[[0, 1,2,3, 4]], axis=1)
    # row wise mean of samples of each metabolite in PC
    PC['mean']=PC1.mean(axis=1)

    LPC1=LPC.drop(LPC.columns[[0, 1,2,3, 4]], axis=1)

    LPC['mean']=LPC1.mean(axis=1)

    Plasmalogen1=Plasmalogen.drop(Plasmalogen.columns[[0, 1,2,3, 4]], axis=1)

    Plasmalogen['mean']=Plasmalogen1.mean(axis=1)
    #Similarity b/w two datasets: MWW ranksum test
    z_stat_1, p_val_pc_lpc = stats.ranksums(PC['mean'],LPC['mean'])      

    z_stat_2, p_val_plas_pc = stats.ranksums(Plasmalogen['mean'],PC['mean'])      

    z_stat_3, p_val_lpc_plas = stats.ranksums(LPC['mean'],Plasmalogen['mean'])      
    #Similarity b/w all the three datasets : ANOVA
    f_val, p_val_oneway = stats.f_oneway(PC['mean'],LPC['mean'],Plasmalogen['mean'])  
    output = {
        "filename":uploadedFile,
        "df":{
            "rows":len(df),
            "columns":len(df.columns),
            "list_of_metabolites":pd.Series(df['Accepted_Compound_ID']).to_json()
        },
        "PC":{
            "rows":len(PC),
            "columns":len(PC.columns),
            "list_of_metabolites":pd.Series(PC['Accepted_Compound_ID']).to_json()
        },
        "LPC":{
            "rows":len(LPC),
            "columns":len(LPC.columns),
            "list_of_metabolites":pd.Series(LPC['Accepted_Compound_ID']).to_json()
        },
        "Plasmalogen":{
            "rows":len(Plasmalogen),
            "columns":len(Plasmalogen.columns),
            "list_of_metabolites":pd.Series(Plasmalogen['Accepted_Compound_ID']).to_json()
        },
        "Retention_time":pd.Series(df['Retention_time_(min)']).to_json(),
        "Retention_time_round_off":pd.Series(df['Retention_time_round_off(min)']).to_json(),
        "mww":{
            'pc_lpc':p_val_pc_lpc,
            'plas_pc':p_val_plas_pc,
            'lpc_plas':p_val_lpc_plas
        },
        "oneway":p_val_oneway

    }
    response=jsonify(output) 
    response.headers.add('Access-Control-Allow-Origin', '*')
    return response

@app.route('/getSimilarMetabolites')
def similar_group(): 
    global uploadedFile   
    global df
    global PC
    global LPC
    global Plasmalogen 
    if(not uploadedFile):
       return redirect(url_for('index'))       
    similar_PC = get_similar_metabolites(PC)
    similar_LPC= get_similar_metabolites(LPC)
    similar_Plasmalogen = get_similar_metabolites(Plasmalogen)    
            
    output = {
        'similar_PC':similar_PC,
        'similar_LPC':similar_LPC,
        'similar_Plasmalogen':similar_Plasmalogen
    } 
    response=jsonify(output) 
    response.headers.add('Access-Control-Allow-Origin', '*')
    return response    

if __name__ == '__main__':
    app.run(
        host="0.0.0.0",
        port=int("8000"),
        debug=True
    )