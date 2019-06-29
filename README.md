# SHAPE_library
Simulated Haploid Asexual Population Evolution (SHAPE)

This is SHAPE modified to be an R library rather than an assemblage of scripts.  I've re-tooled the script set with which I developed SHAPE (found here: {https://github.com/JDench/SHAPE}) so that they form this R-library/package with the improved runtime offered by compilation.  There are only three functions that an average user would really ever need, these are:

defineSHAPE() - Call this after loading the library as it will initialise a set of parameters (as global options).  Once familiar with SHAPE, and ready to manually change parameters, you can call this function while passing arguments in order to override defaults.
designExperiment() - Call this if you want to build an in silico experiment for a breadth of parameters.  Read more about the function to understand the arguments to pass but, in short you'll find in this Git repository there are the SHAPE_templates, in which you'll find the experimentalDesign.v.#.r file where you can define parameters which will be used to create your experiment.  This experimental design will create directories on your machine as well as shell scripts to automatically run the experiment.  Enjoy!
runSHAPE() - This is the actual run function and is a wrapper to call almost everything else defined in the library.  If you just want to see that SHAPE works, load the library, run defineSHAPE(), then run this and you'll get a default experiment run on your machine with output found in a new directory located under your current relative filepath.  

I hope you'll find SHAPE useful for supporting theoretical and empirical study of evolution in an asexual system. I'll note that I've done my best to make SHAPE robust to errors but it's not perfect. If you intentionally try to input parameters that will break the program I have no doubt you can.  However, if you're intent is to use the tool to study meaningful evolutionary parameters, including genomes with millions of sites, communities with population sizes in the billions, and a breadth of fitness landscape models, well - I reckon you could be very pleased.

Again, thanks for choosing to use SHAPE, please feel free to contact me if you have any questions or concerns:

Sincerely,

Jonathan Dench Jdenc017@gmail.com

############################################################################## ############################## CITATION ###################################### ############################################################################## If you plan to cite SHAPE the related manuscript is pending peer review submission but can be found on BioRxiv:

Dench, Jonathan (2018) The SHAPE of logistic growth shows that timing does matter, bioRxiv, doi:10.1101/392944

BibTex Citation Script:

@article {Dench:2018aa, author = {Dench, Jonathan}, title = {The SHAPE of logistic growth shows that timing does matter}, year = {2018}, doi = {10.1101/392944}, publisher = {Cold Spring Harbor Laboratory}, abstract = {Experimental evolution is a powerful tool for studying evolutionary questions and the use of in silico systems is growing as they offer complete parameter control and more fine grained tracking of dynamics. However, existing software are limited by the models implemented, output obtainable, or their lack of interpretability within biological systems. Here I introduce SHAPE (Simulated Haploid Asexual Population Evolution) a forward time simulation tool closely replicating microbial experimental evolution that offers the flexibility to implement alternative models and study detailed output. I first validate SHAPE by replicating the evolutionary outcome expected by several theoretical works. I then extend theory by contrasting how serial passaging affects the evolutionary outcome for microbial communities undergoing exponential and logistic growth.}, URL = {https://www.biorxiv.org/content/early/2018/08/16/392944}, eprint = {https://www.biorxiv.org/content/early/2018/08/16/392944.full.pdf}, journal = {bioRxiv} }
