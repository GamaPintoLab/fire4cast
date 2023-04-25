# fire4cast
Rule Based Model to predict risk of Fire blight disease transmission 

This repositories contains the files necessary to:

1 - Reproduce the training of model rules to predict Fire blight disease transmission based on data collected in two pear orchards in the Alcobaça region (Portugal) during the period February-June between 2019 and 2022.  

Run the R script "f4c_script_2git.R" for which you will need the files with climate data ("vn_qn_clim.xlsx"), phenotype state data ("feno2r.xlsx"), bacterial detection data ("sampling2r.xlsx") and the information at which sampling dates an increase in infection was observed ("alarms.txt").

2 - Use the trained rule set to make new predictions 

Run the R script "fire4cast_new_pred_script.R" for which you will need the file with the model rules ("rulemat4.csv") and the new climate data (an example is "clim_example.csv"). You should prepare a new climate data with the same format as "clim_example.csv". Each climate variable should be the accumulated sum of the 4 days before the prediction date.
