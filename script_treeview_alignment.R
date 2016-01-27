library(ggtree)
fasta <- ("process_AtUGTPROT.fasta.final_tree.fa")
best_tree <- read.best("process_AtUGTPROT.fasta.final_tree.nw")
msaplot(ggtree(tree, size = .5, ledderiza=FALSE)+geom_text(aes(label=label), size=3, color="purple", hjust= 0 ), fasta)



tree <- read.tree("process_AtUGTPROT.fasta.final_tree.nw" )
ggtree(tree)

fasta <- system.file("examples/FluA_H3_AA.fas", package="ggtree")
msaplot(ggtree(tree), fasta) 