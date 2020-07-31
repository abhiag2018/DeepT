import preptools as pt
base_dataDir = "/projects/li-lab/agarwa/CUBE/DeepTact/dataset"
EP_dataPath = 'Javierre_ref_18/DATA_S1/PCHiC_peak_matrix_cutoff5.tsv'
for cell in ['tB','tCD4','nCD4','FoeT','Mon','tCD8']:
	PE_data = pt.training_PE(base_dataDir,EP_dataPath,selectCell=cell,threshold=5,intTypeList=['PE','PP'],output_name='PCHiC')
