library( sarscov2 ) 
library( ape ) 
library( lubridate )
library( treedater )

#' For a given COGUK region: Estimate time trees, compute spike 614 genotype, compute coalescent estimates of growth rates over time for spike D & G
#' 
#' Note this requires development versions of skygrowth and treedater 
#' 
#' @param r The region corresponding to an element in the majora 'adm2' field. Can also be a regular expression to match multiple elements in adm2. 
#' @param tr1 A maximum likelihood phylogeny. Can contain more tips than are in the region 'r', in which case these will be pruned out 
#' @param maj A data frame with the majora metadata. This **must** contain an extra column 'tip.label' which matches labels in given ML tree 
#' @param cogfasta A path to a fasta alignment. This **must** have names that correspond to tips in ML tree 
#' @param ncpu Number of CPUs 
#' @param ntres Number of distinct time-scaled phylogenies to estimate based on ML tree. Coalescent inference will be averaged over these. 
#' @param res Number of coalescent growth rate change points 
#' @param tau0 Scale of exponential prior on change in growth rate 
#' @param ... Additional args are passed to skygrowth1 ( e.g. tstart )
#' @return A list with elements:
#' - spike: Contains coalescent estimates of Ne and growth for D and G 
#' - tds: List of treedater time trees 
#' - data: sub-alignment 
#' - s614: Data frame with spike 614 genotypes 
spike0 <- function(r='BUCKINGHAMSHIRE'
  , tr1 = NULL 
  , maj = NULL 
  , cogfasta = NULL
  , ncpu = 6
  , ntres = 10 
  , res = 10 
  , tau0 = 30 / 365/ 36.5^2
  , ... 
){
# check input: 
stopifnot( 'tip.label' %in% colnames(maj )  )
stopifnot( all( tr1$tip.label %in% maj$tip.label  )  )
	
	i <- which( grepl( maj$adm2, patt = r ) )
	maji <- maj[i, ]
	
	maji$sts <- decimal_date( as.Date( maji$collection_date ) )
	maji <- maji[ !is.na(maji$sts) , ]
	
	tr2 = keep.tip( tr1, maji$tip.label )
	tr2 = unroot(multi2di(tr2))
	
	d = read.dna( cogfasta , format = 'fasta' )
	# check aligment names 
	stopifnot( all( tr2$tip.label %in% rownames(d)   )  )
	
	# make time trees 
	maji <- maji[ match( tr2$tip.label, maji$tip ) , ]
	maji$tip.label2 = paste( sep = '|' , maji$tip.label, maji$sts, '_Il' )
	tr2$tip.label <- maji$tip.label2
	tds = make_starting_trees( fastafn = tr2, ncpu = ncpu , ntres = ntres , treeoutfn = paste0( r, '.nwk') ) 
	# smooth node times 
	gtds = parallel::mclapply( tds, function(td) gibbs_jitter( td, returnTrees=2 )[[2]] , mc.cores = ncpu )
	
	# make sub alignment 
	di = d[maji$tip.label ,] 
	rownames(di ) = maji$tip.label2[ match( rownames(di), maji$tip.label ) ]
	algnfn = tempfile() 
	write.dna( di, file = algnfn, format = 'fasta' )
	
	# compute genotype 
	s614 = compute_spike614genotype(  algnfn )
	
	# do the thing 
	sp = s614_phylodynamics( gtds, s614 , res = res, tau0 = tau0, ...  )  
	
	list( spike = sp , tds = gtds, data = di , s614 = s614 , region = r )
}


# Example: 
if (F)
{
	# For a quick start, run this code to use pre-computed ML tree 
	tr0 = read.tree( 'cog_global_2020-05-15_tree.newick' ) 
	
	# Load the metadata, jump through some hoops to match with ML tree labels 
	majora = read.table( 'majora.20200515.metadata.tsv', header=TRUE, stringsAs=FALSE , sep = '\t')
	majora$id3 = sapply( strsplit( majora$secondary_identifier , 'hCoV-19/'), function(x) x[2] ) 
	majora$id4 = sapply( strsplit( majora$id3 , '/'), function(x) x[2] ) 

	itip2maj <- match( sapply(strsplit( tr0$tip, '/'), '[', 2), majora$central_sample_id )
	itip2maj2 = match( sapply(strsplit( tr0$tip, '/'), '[', 2), majora$id4 )
	itip2maj[is.na(itip2maj)] <- itip2maj2[ is.na( itip2maj) ]
	
	# put this in the same order as the tree tip labels 
	maj = majora[ itip2maj, ]
	
	# REMEMBER to define the tip.label field 
	maj$tip.label <- tr0$tip.label 
	
	# drop anything that couldn't be matched 
	tr1 <- drop.tip( tr0, tr0$tip[ is.na( itip2maj ) ] )
	maj <- maj[ match( tr1$tip , maj$tip ) , ]
	
	# Run the analysis 
	o = spike0( maj = maj
	 , tr1 = tr1 
	 , cogfasta = 'cog_2020-05-15_all_alignment.fasta'
	 , r = 'EDINBURGH'
	) 
	
	# plots
	pne = with( o$spike, add_ggNe_sarscov2skygrowth( sgD, sgG ) )
	pgr = with( o$spike, add_ggGR_sarscov2skygrowth( sgD, sgG ) )
	spl = plot_sample_distribution_spike614( o$s614 )
	
	print( pne )

}

