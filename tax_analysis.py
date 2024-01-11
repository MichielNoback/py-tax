import pandas as pd
from model import Tree

def load_count_data(counts_file):
    counts = pd.read_table(counts_file, index_col=0)
    return counts

def normalize_counts(counts_df):
    counts_columns = counts_df.columns.difference(['OTUID'])
    normalized_counts = counts_df.copy()

    for count_column in counts_columns:
        normalized_counts[count_column] = normalized_counts[count_column] / normalized_counts[count_column].sum() * 100

    return normalized_counts

def process_counts(tree, normalized_counts, cutoff_percentage=1):
    for zotu in normalized_counts.index:
        for sample in normalized_counts:
            # if sample != 'S001P8292':
            #     continue
            count = normalized_counts.loc[zotu, sample]
            if count < cutoff_percentage:
                continue
            else:
                node = tree.zotus[zotu]
                node.counts[sample] = count
                parent = node.parent
                #print(f'zotu:{zotu}; sample={sample}; name:{node.name}; zotus:{node.zotus}; parent:{parent.name}; count:{round(node.counts[sample], 2)}')
                while parent.name != 'root':
                    node = parent
                    parent = node.parent
                    node.counts.setdefault(sample, 0)
                    node.counts[sample] = node.counts[sample] + count
                    #print(f'\tname:{node.name}; zotus:{node.zotus}; count:{round(node.counts[sample], 2)}; parent:{parent.name}')

def get_counts(sample, tax_level):
    print(f'counts for {sample} at level {tax_level}')
    #sum = 0
    counts = {}
    for node in tree.nodes.values():
        if node.level == tax_level and sample in node.counts:
            #sum += node.counts[sample]
            counts[node.name] = node.counts[sample]
            #print(f'node {node.name} has count {node.counts[sample]}')
    #print(sum)
    return counts

if __name__ == "__main__":
    tree = Tree()
    tree.add_lineages_file('ASV_taxonomy_extra_small.txt')
    #print(str(tree))
    # for zotu in tree.zotus:
    #     print(f'zotu:{zotu}; name:{tree.zotus[zotu].name}; zotus:{tree.zotus[zotu].zotus}')

    counts_file = "ASV_counts_extra_small.tsv"
    counts_df = load_count_data(counts_file)
    normalized_counts = normalize_counts(counts_df)
    #print(normalized_counts)

    process_counts(tree, normalized_counts)

    counts = get_counts('S002P8292', 'phylum')
    print(counts)
