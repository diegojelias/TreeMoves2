import dendropy

t1 = "((((((((((((((((((((((((a,b),c),d),e),f),g),h),i),j),k),l),m),n),o),p),q),r),s),t),u),v),w),x),y);"



t2 = "(((((((((((((((((((((((b,c),d),e),f),g),h),i),j),k),l),m),n),o),p),q),r),s),t),u),v),w),x),(a,y));"



taxa = dendropy.TaxonNamespace() #set taxa same for all 
tree1=dendropy.Tree.get(data=t1,schema='newick',taxon_namespace=taxa, rooting='force-rooted')
tree2=dendropy.Tree.get(data=t2,schema='newick',taxon_namespace=taxa, rooting='force-rooted')
tree1.encode_bipartitions()
tree2.encode_bipartitions()
dist=dendropy.calculate.treecompare.symmetric_difference(tree1,tree2)
dist