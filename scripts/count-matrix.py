import pandas as pd

counts = []
for f in snakemake.input:
    coverage = pd.read_table(f, header=None, names=['chr', 'start', 'end', 'coverage'])
    coverage['chr'] = coverage['chr'].astype(str) + "." + coverage['start'].astype(str) + "." + coverage['end'].astype(str)
    coverage.set_index(['chr'], inplace=True)
    counts.append(coverage[['coverage']])

for t, (sample) in zip(counts, snakemake.params.samples.index):
    t.columns = [sample]

matrix = pd.concat(counts, axis=1)
matrix.index.name = "peak"
# collapse technical replicates
#matrix = matrix.groupby(matrix.columns, axis=1).sum()
print(matrix)
matrix.to_csv(snakemake.output[0], sep="\t")
