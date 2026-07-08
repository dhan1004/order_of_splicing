
from os.path import join
import xml.etree.ElementTree as xml_tree
import re
from collections import Counter
import h5py
import gzip

# This scripts filters the GEO metadata associated with high-throughput sequencing studies in Homo sapiens
# in order to find studies that fit certain criteria (e.g. study design, cell line was used, etc.) and
# don't mention sample characteristics we want to exclude (e.g. polyA selection, single cell RNA-seq) 

if __name__ == '__main__':

    project_dir = '/users/dhan30/splicing_order'
    data_dir = join(project_dir, 'data')
    gse_xml_hdf5_path = join(data_dir, 'mouse_expression_profiling_by_HT_seq_gds_metadata.hdf5')
    gse_xml_hdf5_indices_path = join(data_dir, 'mouse_expression_profiling_by_HT_seq_gds_metadata_indices_3t3.txt')
    gse_candidates_file = join(data_dir, 'mouse_GSE_candidates_3t3.txt.gz')

    # I was looking for hek293 total RNA-seq so I used regular expressions to find these terms in the Series Design metadata
    # You can alter these to look for other information, e.g. specific cell lines
    # That sort of information may also be present under the Sample metadata within the Characteristics metadata field
    # so you may want to create a regular expression that searches there
    experiment_include_keywords = [r'3t3']

    # You should keep these as single cell data will not be suitable for the analysis you're doing
    experiment_exclude_keywords = [r'scRNA', r'single-cell', r'single cell', r'sc-RNA']
    experiment_include_re = r'|'.join(experiment_include_keywords)
    experiment_exclude_re = r'|'.join(experiment_exclude_keywords)
    
    # This regular expression should ideally remove polyA-selected samples by using the info in the Extraction Protocol metadata
    # It's not perfect though as these metadata is sometimes incomplete or incorrect
    extraction_exclude_keywords = [
        r'poly\WA', r'polyA', r'poly-A', 
        r'oligo\WdT', r'oligo-dT', r'oligo\(dT\)',
        r'mRNA-seq', r'mRNA seq', r'mRNA librar',
        r'messenger RNA', r'mrna',
        r'poly\(A\)', r'polyA\+',
        r'3\' bias', r'3\'-bias',
        r'3\' end', r'3\'-end'
    ]
    extraction_exclude_re = r'|'.join(extraction_exclude_keywords)

    # Two files will be output:
    # 1. A file of candidate GSE studies that pass the regular expression filters defined above 
    #    as well as some other metadata filters (e.g. "extracted_molecule == 'total RNA' and lib_strategy == 'RNA-Seq' and lib_source == 'transcriptomic'")
    # 2. An index file which will associate the GSE study IDs with there index in the hdf5 file. The hdf5 file is essentially just a large array
    #    that associates an index with a string containing the entire XML metadata file for that GSE study.
    with h5py.File(gse_xml_hdf5_path, 'r') as in_file, gzip.open(gse_candidates_file, 'wt') as out_file, open(gse_xml_hdf5_indices_path, 'w') as index_out:
        out_file.write(f'GSE_ID\tmatch\texample_extraction_description\n')
        index_out.write('GSE_ID\thdf5_index\n')
        
        # This is a prefix present in the titles of metadata elements for some reason
        # so it is needed when searching for the elements in the XML tree
        xmlns = '{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}' 
        series_counter = Counter()
        for i, xml_string in enumerate(in_file['GSE_XML_strings']):
            if len(xml_string) == 0:
                continue
            series_counter['series_processed'] += 1
            
            # Convert the XML string from the hdf5 to a tree structure
            gse_xml_tree = xml_tree.fromstring(xml_string)

            # Get the series node and write the GSE ID index to the index file
            series_node = gse_xml_tree.find(f'{xmlns}Series')
            gse_id = series_node.attrib['iid']
            index_out.write(f'{gse_id}\t{i}\n')
            series_title = series_node.find(f'{xmlns}Title').text.strip()
            series_summary = series_node.find(f'{xmlns}Summary').text.strip()
            series_design = series_node.find(f'{xmlns}Overall-Design').text.strip()

            # Only proceed if there are no matches in the Overall Design metadata to the regular expression defined by experiment_exclude_re
            if re.search(experiment_exclude_re, series_design, flags=re.IGNORECASE) is None:
                experiment_list = []
                # This searches for the regular expression defined by experiment_include_re
                # You can then examine the match to pull out relevant information, e.g.
                # if you searched for a list of cell lines, which cell line matched?
                match_list = []
                for experiment_match in re.finditer(experiment_include_re, series_design, flags=re.IGNORECASE):
                    match_start, match_end = experiment_match.start(), experiment_match.end()
                    match_str = series_design[match_start:match_end]
                    match_list.append(match_str)
                
                if len(match_list) > 0:
                    matches_combined = ','.join([m.replace('\n', ' ') for m in match_list])         
                    series_samples = gse_xml_tree.findall(f'{xmlns}Sample')
                    sample_count = len(series_samples)

                    extraction_pass, lib_pass = True, False
                    paired_end_samples = 0
                    # The sample node is where you would find the Characteristics metadata which might also contain info on cell lines
                    for sample_node in series_samples:
                        lib_strategy = sample_node.find(f'{xmlns}Library-Strategy').text.strip()
                        lib_source = sample_node.find(f'{xmlns}Library-Source').text.strip()
                        layout_node = sample_node.find(f'{xmlns}Library-Layout')
                        if layout_node is not None:
                            is_paired = layout_node.find(f'{xmlns}PAIRED') is not None
                            if is_paired:
                                paired_end_samples += 1
                         # Check read length in characteristics
                        characteristics_nodes = sample_node.findall(f'{xmlns}Characteristics')
                        for char_node in characteristics_nodes:
                            if char_node.text:
                                char_text = char_node.text.strip().lower()
                                if 'read length' in char_text or ('bp' in char_text and 'read' in char_text):
                                    length_match = re.search(r'(\d+)\s*bp', char_text)
                                    if length_match and int(length_match.group(1)) < 80:
                                        extraction_pass = False
                                        break

                        for channel_node in sample_node.findall(f'{xmlns}Channel'):
                            extracted_molecule = channel_node.find(f'{xmlns}Molecule').text.strip()
                            extraction_protocol = channel_node.find(f'{xmlns}Extract-Protocol').text.strip()
                            # Check that there are no matches to the regular expression meant to filter out mentions of polyA
                            if re.search(extraction_exclude_re, extraction_protocol, flags=re.IGNORECASE) is not None:
                                extraction_pass = False
                            if extracted_molecule == 'total RNA' and lib_strategy == 'RNA-Seq' and lib_source == 'transcriptomic':
                                lib_pass = True

                    # require majority paired end samples
                    paired_end_ratio = paired_end_samples / sample_count
                    if paired_end_ratio > 0.75:  # At least 80% paired-end
                        continue
                        
                    if extraction_pass and lib_pass:
                        # Output the GSE IDs for studies that pass all filters along with their match strings and an example extraction protocol string
                        # You could manually examine the extraction protocol to confirm it doesn't appear to polyA/mRNA
                        # out_file.write(f'{gse_id}\t{matches_combined}\t{extraction_protocol.replace("\n", "\t")}\n')
                        extraction_protocol_clean = extraction_protocol.replace("\n", "\t")
                        out_file.write(f'{gse_id}\t{matches_combined}\t{sample_count}\t{paired_end_ratio:.2f}\t{extraction_protocol_clean}\n')
                        series_counter['candidates_found'] += 1
                        if series_counter['candidates_found'] % 10 == 0:
                            print(f'{series_counter["candidates_found"]} candidate GEO Series found out of {series_counter["series_processed"]} series...')
                    
                print(f'\nFinal results: {series_counter["candidates_found"]} stringent candidates from {series_counter["series_processed"]} total series')
                print(f'Reduction factor: {series_counter["series_processed"] / max(series_counter["candidates_found"], 1):.1f}x')
                




