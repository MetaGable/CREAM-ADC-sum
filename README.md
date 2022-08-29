# CREAM L1 data

### Goals
Our goal is to produce L1 data (convert ADC into Energy units). We are currently deriving the gain for mid ranged channels.   

### Notes on our scripts 
######## Most if not all requires extensive adjustment to execute in desired directory on desired data.
######## Most sample of the plots generated could be found in my presentation on the ELMS page, though the display style may vary as the scripts updates. 
###### ScatterV8.py
 - Create data files for Plot_tackedV2.py to plot. It organizes event signal based on ribbons and creates pairs containing low range and mid range ADC. 
 - This is the first step of the analysis and will need to be modified and executed again when analyzing mid & high range ADC. 
 
###### Plot_all_evt.py
 - The second step. Plot and fit the result from scatter for both full and recontructed (aka tracked) data. 
 - Like most of the plots, you'll comment in and out to change the type of data it is feeding on. 
 - One will need to adjust some parameters in the script, particularly x-factor and y-factor according to the data. different data type may require different parameters. 

###### manual.py
 - The thrid step. For manually specifing the cuts when automatic cut in Plot_all_evt.py fails. 
 - Refer to full_manual_record2.csv for the format of specifing cuts. 
 
###### time_gradient.py
 - Jointly the thrid step. Used to roughly identify when was the event recorded. 
 - For mac uses, there is an app called Digital Color Meter is used to help interpret the plot. 
 - Refer to full_manual_record_cloud.csv for the format of specifing cuts. 
 - Unfortunately I don't have a good script for removing data from some time period in a scatterplot. One will have to modify this code to do so. 
  
###### compile_table.py
 - The final step, used to gather your results into a final table when use finished generating all the plots and made sure all the fits are good.  

#### Outdated and kept for documentation purposes
###### Plot_tackedV2.py
 - Plot and fit the result from scatter for recontructed data. 

###### TG.csv
 - The outdated input format for time_gradient.py

###### manual_dg.py
 - Not exactly outdated but an alternate version of the manual.py code above.
 
