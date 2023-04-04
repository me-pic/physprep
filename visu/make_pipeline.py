"""
Movation
--------
Create a function that will generate a figure containing the processing steps with the type of filter used, 
and the value of those parameters. This function will enable users to clearly show their processing pipeline
in order to increase the reproducibility of their signal processing. 
"""
def make_pipeline(pipeline_dict):
"""
One way to do it:
if pipeline = True in neuromod_process.py function, create a dictionnary that will be incremented after each processing step
The key will be the step number (ex. 1, 2, 3, 4, ...), and the value will be another dictionnary. That dictionnary will contain
the filter name (key 'filter'), the filter parameter (key 'parameters'), and a 10s window of the signal after that filter applied 
(key 'signal')
pipeline_dict = {
    '1': {
        'Filter': str,
        'Parameters': {}
        'Signal': ndarray
    },
    '2': {

    },
    ...
    'n': {

    }
}
"""
