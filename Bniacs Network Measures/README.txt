README.txt

BNIACS UPDATED FOLDER:

These are explanations for what each .py file.


parseNewData.py - This parses a ms.pkl, which is multiple sclerosis data.
    
    The ms.pkl's data structure is a dictionary.
    Key: a two-element tuple. 
        First element referrs to the type of matrix that we have as the value.
        The first element is one of the following: 
              ['corr','lcorr,'lacorr','pcorr','pacorr']
        Note: the ms.pkl file we currently have also has
              ['tcorr','tlcorr','tlacorr','tpcorr','tpacorr']
        but does not currently parse for these values.
              
        Second element corresponds to the groups that patients are grouped in.
        In this case, there are four possible groups:
              ["PP","CIS","RR","SP"]
              
    Value: a list of matrices
        for the MS data, the matrices are of size: 70x70
        
        data = load(open("ms.pkl","rb"))
        data[('corr','PP')] returns a list of matrices corresponding to the
                            patient group: 'PP' and corr matrices.
        
        data[('tcorr','PP')] returns a list of matrices that are _thresholded_.
        _thresholded_ means if we do:
            g = nx.DiGraph(data[('tcorr','PP')][0])
            # g will be a directed unweighted graph.
    
    Currently, this parses the current groups we have.
    If we have more types of matrices, add to the list in line 12.
    If we have more group of patients, add to the list in line 13.
    Once you have an updated dictionary, change the FNAME in line 14 to the filename of your dictionary.
    
    main1() computes all the raw measures. (NOTE: Creates a bunch of .pkl files in same directory)
    main2() computes all the t-test and creates T-test charts in the same directory.

parseHelper.py - helper functions for parseNewData.py    
    
nalz_test.py - This makes the t-test charts work. This .py file is used by parseNewData.

    If you want to use this just to make t-test charts, use the: mycompare function.
    
    mycompare - Makes the t-test charts and saves them in same directory.
    
        T-test chart comprises of: Each group compared with one another + each group compared with random-stats.
    
        First arg: is a list of stat measures for each group.
    
            for example, I included in the ms_stats folder a bunch of convertedFormat .pkl files.
            
            d = load(open("ms_stats\ms_convertedFormat_pacorr_undirected.pkl","rb"))
            
            len(d) in this case is 4. (There are 4 patient group types)
            
            Each element in d is of the following data structure:
            [('StronglyConnected','avgoutdeg','stdoutdeg, .... ),([1.0,1.0,1.0,....],[10.0,10.0,...] ...)]
            
            ^ Basically a list of two tuples. 
                First tuple is comprised of all the measure-strings.
                Second tuple is comprised of a list of measure-values correponding to the meausre-strings (index correspondence)
                
            Look in ms_stats\ directory for sample .pkl files.
            
        Second arg: a list of the strings that correspond to the groups in the first argument.
        NOTE: must be in the same order. For example if first arg is 
                [stats for 'PP',stats for 'CIS',stats for 'RR',stats for 'SP']
                Then the second arg must be:
                ['PP','CIS','RR','SP']
        
        
        Third arg: File for randomized data (randomly generated matrices with stats computed).
            This argument is the same format as an element in d (line 52 of this README).
        
        (look in ms_stats\uprand_D10.pkl and ms_stats\uprand_D10.pkl)
        
        
        
        Fourth arg: The name that you want to name the t-test .png file. 
        
        NOTE: this will save in the same directory, so if you have a file with same name, WILL OVERWRITE AND REPLACE.
        
mynetalgs.py - This file computes MEASURES.
    
    myallmeasures - computes measures of a simple directed (or undirected) graph.
        input: directed or undirected graph.
        
        output: a tuple. 
            First element: ['stronglyConnected','avgoutdeg',...]
            Second element: [.83,1.0,0.0,...]
            
            Both lists are equal in length.
        
greedy.py - Helper functions for mynetalgs' calculation of modularity and linkrank.

sepGroups.py - seperates a PATIENT ID dictionary into the 4 groups NL,AD,MCI, and CONVERTED.

    Structure of Patient ID dictionary:
            Key: PATIENT ID (string)
            
            Value: Three element tuple:
                First element: Dictionary of measures 
                    key: 2 element tuple:
                            first element: String of the measure.
                            second element: String of the type of corr.
                    value: a number.
                
                Second element: Dictionary of matrices
                    key: one of the following - ['corr','lcorr,'lacorr']
                                    or ['tcorr','tlcorr','tlacorr']
                                (t means thresholded)
                    
                    value: a matrix.
                    
                Third element: A string that represents the patient's group.
                    Either "AD" "MCI" "NL" or "CONV".
    
    To seperate this patient ID dictionary into groups: call the sepGroups function with the string
        of your patient ID file name.
        
    Will save as .pkl files in same directory named adCorr,mciCorr,nlCorr,convCorr.
        ^ are basically a list of matrices.
        
corrData.py - fills in the measures for a PATIENT ID dictionary with an empty dictionary of measures.

    main() is the function that does this.
    
    First argument: string of the file name of the PATIENT ID dictionary .pkl file (without the measures)
    
    Second argument: string of the file name that you want the new measure-filled dictionary to me written as.
    
    For example:
    main("aznorbert_corrsd_undirected.pkl","aznorbert_corrsd_new_measures_undirected.pkl")
    
DataTest.py - tests whether two the corrData.py worked. Compares two PATIENT ID dictionaries.
    
    -replace line 10 and 14 with appropiate file names.
    
graph_helper.py - function for graph-swapping. Not used anywhere else.
    

    
    
    
                    
                    
                    
                
        
                    
            



    
    
    
    
    




