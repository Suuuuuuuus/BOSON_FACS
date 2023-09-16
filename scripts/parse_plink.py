import FACSus as fs

markers = pd.read_csv('markers.txt', header = None).rename(columns = {0: 'marker'})
markers = markers.marker.to_list()
for i in markers:
    fs.parse_plink('results/' + i + '.assoc.linear', True, 'results/' + i + '.baseline.gwas.txt')
