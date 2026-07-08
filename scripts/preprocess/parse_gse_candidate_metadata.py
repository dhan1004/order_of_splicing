from os.path import join
import xml.etree.ElementTree as xml_tree
from collections import Counter
import h5py
import gzip

# This script takes the candidate GSE studies output by filter_gse_miniml_metadata.py and 
# pulls the SRA run IDs associated with each of the samples in the study. It outputs a new
# hdf5 file that contains metadata strings in a simpler format than those in the original
# GEO XML file.

if __name__ == '__main__':

	project_dir = '/users/dhan30/splicing_order'
	data_dir = join(project_dir, 'data')
	gse_xml_hdf5_path = join(data_dir, 'mouse_expression_profiling_by_HT_seq_gds_metadata.hdf5')
	gse_xml_hdf5_indices_path = join(data_dir, 'mouse_expression_profiling_by_HT_seq_gds_metadata_indices_3t3.txt')
	gse_candidates_file = join(data_dir, 'mouse_GSE_candidates_3t3.txt.gz')
	gse_candidates_sample_metadata_file = join(data_dir, f'mouse_gene_perturbation_GSE_candidates_sample_metadata_3t3.hdf5')
	gse_candidates_sample_info_file = join(data_dir, f'mouse_gene_perturbation_GSE_candidates_sample_title_characteristics_3t3.txt.gz')
	gsm_sra_pipeline_file = join(data_dir, 'mouse_gsm_sra_list_for_pipeline_3t3.tsv')
	sra_gsm_table = join(data_dir, 'SRA_GSM_accessions.txt.gz')

	# Parse the association between GSE ID and the index of its metadata string in the original hdf5 file
	gse_indices = {}
	with open(gse_xml_hdf5_indices_path, 'rt') as in_file:
		next(in_file)
		for line in in_file:
			gse_id, index = line.strip().split('\t')
			gse_indices[gse_id] = int(index)
	print('Done parsing GSE hdf5 indices...')
	
	# Parse the association of GSM ID to the ID(s) of its SRR run(s)
	gsm_to_sra = {}
	with gzip.open(sra_gsm_table, 'rt') as in_file:
		for line in in_file:
			srr_id, gsm_id = line.strip().split('\t')
			gsm_id = gsm_id.split('_')[0].split('.')[0]
			if gsm_id not in gsm_to_sra:
				gsm_to_sra[gsm_id] = []
			gsm_to_sra[gsm_id].append(srr_id)
	print('Done parsing GSM to SRR table...')

	# Parse the GSE IDs of the candidate studies
	gse_id_list = []
	with gzip.open(gse_candidates_file, 'rt') as in_file:
		next(in_file)
		for line in in_file:
			gse_id = line.strip().split('\t')[0]
			gse_id_list.append(gse_id)
	
	# Output two files:
	#  1. A file that associates each GSE ID with all sample GSM IDs along with their sample titles and sample characteristics
	#  2. A new hdf5 with custom metadata strings for each GSE ID where each string comprises the GSE ID, series title,
	#     series summary, series design, and then a substring for each sample that contains its GSM ID, sample title,
	#     and SRA IDs (these could either be SRR (run) IDs or SRX (experiment) IDs)
	xmlns = '{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}'
	sample_counts = Counter()
	with h5py.File(gse_xml_hdf5_path, 'r') as hdf5_file, h5py.File(gse_candidates_sample_metadata_file, 'w') as out_file, gzip.open(gse_candidates_sample_info_file, 'wt') as sample_out, open(gsm_sra_pipeline_file, 'wt') as pipeline_out:
		sample_out.write('gse_id\tgsm_id\tsample_title\tsample_characteristics\n')
		pipeline_out.write('gsm_id\tsra_ids\n')
		str_datatype = h5py.string_dtype(encoding='utf-8')
		str_dataset = out_file.create_dataset('GSE_sample_info_strings', shape=(len(gse_id_list),), dtype=str_datatype, compression='gzip')
		for i, gse_id in enumerate(gse_id_list):
			xml_string = hdf5_file['GSE_XML_strings'][gse_indices[gse_id]]
			gse_xml_tree = xml_tree.fromstring(xml_string)
			series_node = gse_xml_tree.find(f'{xmlns}Series')
			series_title = series_node.find(f'{xmlns}Title').text.strip()
			series_summary = series_node.find(f'{xmlns}Summary').text.strip()
			series_design = series_node.find(f'{xmlns}Overall-Design').text.strip()
			series_samples = gse_xml_tree.findall(f'{xmlns}Sample')
			sample_info_strings = []
			for sample_node in series_samples:
				lib_strategy = sample_node.find(f'{xmlns}Library-Strategy').text.strip()
				lib_source = sample_node.find(f'{xmlns}Library-Source').text.strip()
				lib_pass, sample_char = False, []
				for channel_node in sample_node.findall(f'{xmlns}Channel'):
					extracted_molecule = channel_node.find(f'{xmlns}Molecule').text.strip()
					organism = channel_node.find(f'{xmlns}Organism').text.strip()
					if organism == 'Mus musculus' and extracted_molecule == 'total RNA' and lib_strategy == 'RNA-Seq' and lib_source == 'transcriptomic':
						lib_pass = True
						for char_node in channel_node.findall(f'{xmlns}Characteristics'):
							if len(char_node.attrib) > 0:
								char_type, char_value = char_node.attrib['tag'], char_node.text
								sample_char.append(f'{char_type}: {char_value}'.strip().replace('\n', '').replace('\t', '').replace(';', ','))
				if lib_pass:
					sample_type = sample_node.find(f'{xmlns}Type').text.strip()
					if sample_type == 'SRA':
						# If the GSM ID is present in the GSM to SRR table, use that. It is sometimes easier to explicitly
						# know all the SRR run IDs for a GSM sample but not all of them appear in this table for some reason
						sample_gsm_id = sample_node.attrib['iid']
						if sample_gsm_id in gsm_to_sra:
							sample_title = sample_node.find(f'{xmlns}Title').text.strip().replace(' ', '_').replace('\n', '').replace('\t', '')
							sample_srr_list = ','.join(gsm_to_sra[sample_gsm_id])
							sample_info_strings.append(f'{sample_gsm_id}\t{sample_title}\t{sample_srr_list}')
							sample_char = '; '.join(sample_char) if len(sample_char)>0 else 'NA'
							sample_out.write(f'{gse_id}\t{sample_gsm_id}\t{sample_title}\t{sample_char}\n')
							pipeline_out.write(f'{sample_gsm_id}\t{sample_srr_list}\n')
							sample_counts['samples_passed'] += 1
						# If the GSM ID is not in the SRR table, parse the SRX ID from the sample node information
						# Each SRX comprises one or more SRR runs.
						else:
							sample_relations = sample_node.findall(f'{xmlns}Relation')
							if len(sample_relations) > 0:
								relation_types = [e.attrib['type'] for e in sample_relations]
								if 'SRA' in relation_types:
									sample_title = sample_node.find(f'{xmlns}Title').text.strip().replace(' ', '_').replace('\n', '').replace('\t', '')
									sample_srx_id = [e for e in sample_relations if e.attrib['type'] == 'SRA'][0].attrib['target'].split('term=')[-1]
									sample_info_strings.append(f'{sample_gsm_id}\t{sample_title}\t{sample_srx_id}')
									sample_char = '; '.join(sample_char) if len(sample_char)>0 else 'NA'
									sample_out.write(f'{gse_id}\t{sample_gsm_id}\t{sample_title}\t{sample_char}\n')
									pipeline_out.write(f'{sample_gsm_id}\t{sample_srx_id}\n')
									sample_counts['samples_passed'] += 1
								else:
									sample_counts['sra_id_not_found'] += 1
							else:
								sample_counts['sra_id_not_found'] += 1
						
				sample_counts['samples_processed'] += 1
			
			if len(sample_info_strings) > 1:
				# Each string in the new hdf5 file will have the different info fields separated by line breaks, so we replace any line breaks
				# that occur in the metadata strings with tabs
				gse_info_string = '\n'.join([gse_id, series_title.replace('\n', '\t'), series_summary.replace('\n', '\t'), series_design.replace('\n', '\t')] + sample_info_strings)
				str_dataset[i] = gse_info_string
				sample_counts['study_passed'] += 1
			else:
				sample_counts['study_failed'] += 1
	print(sample_counts)