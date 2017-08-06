# We'll render HTML templates and access data sent by POST
# using the request object from flask. Redirect and url_for
# will be used to redirect the user once the upload is done
# and send_from_directory will help us to send/show on the
# browser the file that the user just uploaded
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

# This route will show a form to perform an AJAX request
# jQuery is loaded to execute the request and update the
# value of the operation
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
        # will basicaly show on the browser the uploaded file
        return redirect(url_for('uploaded_file'))

# This route is expecting a parameter containing the name
# of a file. Then it will locate that file on the upload
# directory and show it on the browser, so if the user uploads
# an image, that image is going to be show after the upload
@app.route('/output')
def uploaded_file():
    return render_template('result.html')

@app.route('/processfile')
def process_file():
    global uploadedFile
    global df 
    global PC
    global LPC
    global Plasmalogen   
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
    global df
    global PC
    global LPC
    global Plasmalogen   
    #create dataframe in such a way that new dataframe is
    #left with sample columns and Accepted_Compound_ID columns
    #which is basically column on metabolites      
    testPC=PC.drop(PC.columns[[0,1]], axis=1)
    #turn the tables i.e
    #Now each entry in col[Accepted_Compound_ID] is now a separate column
    #i.e each metabolite is a separate column now having values as information
    #of 1050 samples
    transPC = testPC.transpose()
    #data cleaning : getting rid of index set as column after transpose
    #setting the first row as column headers
    transPC.columns = transPC.iloc[0]
    transPC = transPC.ix[1:]
    #this can be easily optimized with look_ahead or a backtrack dictionary
    similar_PC=[]
    #iterate thriugh columns of PC,LPC,Plasmalogen
    #to find the MWW score with each metabolite
    for column in transPC.columns:    
        for head in transPC.columns:        
            if column!=head:
                z_stat, p_val = stats.ranksums(transPC[column],transPC[head])
                #if score is greater than 0.05 
                # => metabolites are statistically simillar 
                if p_val>0.05: 
                    arr={
                        "metabolite1":column,
                        "metabolite2":head,
                        "score":p_val
                    } 
                    if len(similar_PC):
                        for item in similar_PC:                        
                            if (item['metabolite1']!=column and item['metabolite2']!=head)and(item['metabolite2']!=column and item['metabolite1']!=head):                                      
                                similar_PC.append(arr)
                    else:
                        similar_PC.append(arr)

    testLPC=LPC.drop(LPC.columns[[0,1]], axis=1)
    transLPC = testLPC.transpose()
    transLPC.columns = transLPC.iloc[0]
    transLPC = transLPC.ix[1:]
    #this can be easily optimized with look_ahead or a backtrack dictionary
    similar_LPC=[]
    for column in transLPC.columns:    
        for head in transLPC.columns:        
            if column!=head:
                z_stat, p_val = stats.ranksums(transLPC[column],transLPC[head])
                if p_val>0.05:
                    arr={
                        "metabolite1":column,
                        "metabolite2":head,
                        "score":p_val
                    }
                    if len(similar_LPC):
                        for item in similar_LPC:                        
                            if (item['metabolite1']!=column and item['metabolite2']!=head)and(item['metabolite2']!=column and item['metabolite1']!=head):                                      
                                similar_LPC.append(arr)
                    else:
                        similar_LPC.append(arr)
    testPlasmalogen=Plasmalogen.drop(Plasmalogen.columns[[0,1]], axis=1)
    transPlasmalogen = testPlasmalogen.transpose()
    transPlasmalogen.columns = transPlasmalogen.iloc[0]
    transPlasmalogen = transPlasmalogen.ix[1:]
    #this can be easily optimized with look_ahead or a backtrack dictionary
    similar_Plasmalogen=[]
    for column in transPlasmalogen.columns:    
        for head in transPlasmalogen.columns:        
            if column!=head:
                z_stat, p_val = stats.ranksums(transPlasmalogen[column],transPlasmalogen[head])
                if p_val>0.05:
                    arr={
                        "metabolite1":column,
                        "metabolite2":head,
                        "score":p_val
                    }
                    if len(similar_Plasmalogen):
                        for item in similar_Plasmalogen:                        
                            if (item['metabolite1']!=column and item['metabolite2']!=head)and(item['metabolite2']!=column and item['metabolite1']!=head):                                      
                                similar_Plasmalogen.append(arr)
                    else:
                        similar_Plasmalogen.append(arr)
            
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
