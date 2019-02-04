## Code to preprocess Detroit Lake empirical data Spring 2019
#! James Watson, The Prediction Lab 2019
import pandas as pd
import numpy as np
import datetime
import openpyxl as px

#### Load excel spreadsheet from Salem
data = px.load_workbook("./PredictionLab_TXN_NUT_YSI_DRAFT.xlsx",data_only=True)
ws = data['PREDICT']

# Extract values
column = ws['E']
for x in np.arange(0,len(column)):
    print(column[x].value)




