# READ ME
This repo contains two main strands of code; a multi-purpose library dataflow to facilitate data management and some biomarkers for (i)EEG seizure severity in the format of the dataflow library. The overview sections below only list the most important files.

The overall rational of this library is that several seizure markers should be calculated quickly and easily. Even when new seizures are added, we should only need to run the analysis on these new seizures. Similarly when new markers are added, their calculation should re-use as much of what is existing already (e.g. if the new marker is based on line-length, you will not need to re-calculate line length again if that was done already for another marker). Thus the following code keeps track of seizures, and calculations (e.g. line length, or bandpower or similar - you can add to it!) in a database, and can use them for subsequent severity measure calculations. 



If you want a quick demo to see what this library can do, start by running **start_here_with_example_data.m**. Some seizure EEG data is included in this repo as a demo, and has been downloaded from the publicly available [ieeg-swez](http://ieeg-swez.ethz.ch/) longterm database.


We are inviting the epilepsy research community to contribute to this library, so please either contact us and/or use pull requests.


19/10/2021, Leonard Waldmann (leonard.waldmann@gmx.de)
July/2022, Sarah Gascoigne (s.gascoigne@newcastle.ac.uk), Yujiang Wang (Yujiang.Wang@ncl.ac.uk)





## Dataflow strand
Dataflow should facilitate working with multiple segments, performing calculations on those segments with different parameters, saving and accessing the calculations by their parametersets, and provide some basic visualisation analysis tools. 

This documentation only summarizes the most important features. To see all functions look into the code. If in need of help, just call help(function_name) or look into the documentation at the beginning of each function in the code.


### Overview
- **Seg_db**: database to store and retrieve segment meta data and/or large data files (e.g. time series)
- **Abstract_calc_db**: abstract class that implements all functionality for computing and storing calculations for different parameter sets. To use functionality, inherit Abstract_calc_db as demonstrated in template_db.m
- **template_db**: template for implementing new calculation database
- **tbl_join**: merge two data tables
- **plot_tbl**: easily plot numeric and/or categorical columns of data tables. Requires the gramm library available on github!
- **gridsearch**: basic parameter gridsearch utility

### Getting started: general pipeline
The general workflow with this library should be: create a calculation database, add some parametersets, calculate the parametersets for some data (eventually retrieved from a Seg_db), get the results from the calculation database, pass them to a ms_ function, join the output from ms_ function with meta data from the data table and visualize results with plot_tbl. Have a look at the pseudo code below and getting_started.m.
```
db = Ampl_db(...)
db.add_paramset(...)
db.calc(seg_tbl,...)
calcs = db.get(...)
ms_tbl = ms_suppr(seg_tbl,calcs,...)
tbl = tbl_join(seg_tbl,ms_tbl)
plot_tbl(...)
```


### Specify segments
Segments have a unique id, the segment_id. Apart from that, they can be specified by their unique index in the database, db_ind, their patient patient_id, or just with the row in the data table. Multiple inputs of the same kind can be combined in numeric-, string-, or cell-arrays. Multiple kinds of inputs can also be combined. If you want to retrieve all data, use a db_ind of -1.
```
seg_db.get_data(‘id’,’GMS21_123’)
seg_db.get_data(‘ind’,4)
seg_db.get_data(‘pat’, {’95’,’1163’})
seg_db.get_data(‘row’,2)
seg_db.get_data(‘row’,20,’pat’,’800’)
seg_db.get_data('ind',-1)
```

### Specify parameter sets
The best way to specify a parameter set (paramset) is a struct that has the parameters as field. After a paramset is added to a database, it also gets a unique integer index in _that_ database that can be used instead of the struct. If a function takes multiple paramsets, their struct or indices can be concatenated in a cell array. If a function takes exactly one paramset, it can either be specified by a positional struct/ind or name value pairs with the parameter names.
```
% multiple parameter sets
ll_db.calc(seg_db,parset_struct)
ll_db.calc(seg_db,3)
ll_db.del_paramset({parset_struct,3})

% single parameter sets
ll_db.add_paramset(’wndw_len’,256,’wndw_overlap’,32)
ll_db.add_paramset(parset_struct)
ll_db.get(2)
```

### General data type
Between functions, data is passed in a table where the columns contain the different data and the rows the different segments. There has to be a column ‘segment_id’ to uniquely identify the rows. Results from different functions can be (inner-) joined with the tbl_join function that matches the rows with identical segment ids.
```
ms = ms_suppr(...);   % calculate measure
ms_tbl = tbl_join(seg_db.get_meta(‘ind’,-1), ms);   % join meta data with measure resuts
```

### Create new calculation database
1. Copy template_db.m, change the class name, change get_algo_params
2. fill in calc_seg: This is the heart of the database and performs the calculation for one segment, if desired include a analysis plot. For detailed info on inputs and required outputs check help(Abstract_calc_db.calc_seg)
3. optional: fill in calc_chn which provides calculation and analysis plots on channel level
Also look at example_calc_db.m.

### Version management
The calculation database also keeps track of the algorithm version as specified in the method get_algo_params. The vers_id is basically treated as additional parameter where vers_cmmt is not compared and just a note for the user. If the requested version id doesn't match the current version id in a calculation request, and error/warning is thrown. However, data from old calculations can further be retrieved by specifying the version id as parameter. So for any major changes in calculation best change the version (or create new database).

### Add new parameters afterwards
After saving calculations in a database, you can still add parameters with a default values for all calculations in the database. Just add them with their default value to get_default_params and call add_defaultparam on your database with the param name and default value (as double check).

### Performance & RAM
When retrieving data from either Seg_db or Abstract_calc_db with the get functions, the data will be loaded from the hard drive into the RAM. If a small subset of data used frequently, it will be way faster to save the retrieved data in a variable in RAM. In Seg_db, this is why there is get_meta and get_data; the first only retrieves the small meta data very fast, where the latter actually loads the whole table from the hard drive. If the ‘heavy’ data is not needed, use get_meta for more speed.

## Biomarker library
This is a collection of potential biomarkers indicating iEEG severity utilizing the dataflow library above. To calculate all measures (suppression measure with amplitude not line length) run the full_analysis script after specifying the input variables.

### Overview
- **Ampl_db**: calculation database (calc db) to store a windowed measure based on amplitude in each channel
- **BP_db**: calc db to store bandpower in each channel
- **DF_db**: calc db to store dominant frequency (DF) timeseries for each channel
- **LL_db**: calc db to store line length (LL)
- **ms_95pctle**: computes the "peak markers" as described in [Gascoigne 2022](https://arxiv.org/abs/2206.15283) that capture the peak amount of a signal feature ever seen in a seizure
- **ms_imprint and ms_propagation**: computes the "imprint markers" as described in [Gascoigne 2022](https://arxiv.org/abs/2206.15283) that capture properties of seizure propagation
- **ms_suppr**: computes post-ictal generalized suppression (PGES) and partial suppression duration and its strength as described in [Gascoigne 2022](https://arxiv.org/abs/2206.15283)
- **ms_dfs_clustering**: [experimental] computes dominant frequency similarity (dfs) for multiple segments, clusters channels by their DFs, and outputs more measures based on the clustering
- **ms_bp_sliding_perio**: [experimental] tries to assign the strength of regular bandpower oszillation in 1.5-2 Hz
- **ms_PLHG**: [experimental] calculated phase locked high gamma properties
