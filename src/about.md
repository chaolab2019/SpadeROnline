---
output: html_document
---

<font color="blue"><h4 align = "right">Original version (March 2015)</h4></font><font color="blue"><h4 align = "right">Latest Version (September 2016)</h4></font>  

<h3 align="center" style = "font-weight: bold">Introduction to Online Program SpadeR</h2>

<h3 align="center" style = "font-weight: bold">(<u>S</u>pecies-richness <u>P</u>rediction <u>A</u>nd <u>D</u>iversity <u>E</u>stimation in <u>R</u>)</h3>

<h5 align="center">by</h5>

<h5 align="center">Anne Chao, K. H. Ma, T. C. Hsieh and Chun-Huo Chiu</h5>

<h5 align="center">Institute of Statistics</h5>

<h5 align="center">National Tsing Hua University, TAIWAN 30043</h5>
================

The program SpadeR is the R-based online version of SPADE available via the <a href="http://chao.stat.nthu.edu.tw/wordpress/software_download/">link</a> or <a href="https://chao.shinyapps.io/SpadeR/">link2</a>. Clicking these links, you will be directed to the online interface window. <font color="red"> Users do not need to learn/understand R to run SpadeR</font>. The interactive web application was built using the Shiny (a web application framework). SpadeR includes nearly all of the important features from the original program SPADE while also having the advantages of expanded output displays and simplified data input formats. Further, owing to the power of the R language, SpadeR now offers high-resolution plots/figures that the original SPADE lacks. Note, however, that some features in SPADE have been expanded to become an independent online program; see below.



Like the original SPADE, the program SpadeR computes various biodiversity indices based on two major types of sample data (abundance data and replicated incidence data) taken from one or multiple communities. A detailed <a href="http://chao.stat.nthu.edu.tw/wordpress/wp-content/uploads/software/SpadeR_UserGuide.pdf">SpadeR User's Guide</a> illustrates how to run this program in an easily accessible way through numerical examples with proper interpretations of portions of the output. SpadeR is divided into six parts:

<li> <u><font color="4464F1">Part I: Species</font></u> (estimating species richness for one community).</li>

<li> <u><font color="4464F1">Part II: Diversity Profile Estimation</font></u>   (estimating a continuous diversity profile in one community including species richness, Shannon diversity and Simpson</br>
&nbsp;&nbsp;&nbsp;diversity). This expanded part also features plots of empirical and estimated continuous diversity profiles. Various estimates for Shannon entropy and the</br>
&nbsp;&nbsp;&nbsp;Gini-Simpson index are also provided.</li>

<li> <u><font color="4464F1"> Part III: Shared Species</font></u>   (estimating the number of shared species between two communities).</li>

<li> <u><font color="4464F1">  Part IV: Two-Community (Similarity) Measures </font></u>    (estimating various similarity indices for two assemblages). The richness-based indices include the </br>
&nbsp;&nbsp;&nbsp;classic two-community Jaccard and Sørensen indices; the abundance-based indices include the Horn, Morisita-Horn, two-community Bray-Curtis and the </br>
&nbsp;&nbsp;&nbsp;abundance-based Jaccard and Sørensen indices.</li>

<li> <u><font color="4464F1">   Part V: Multiple-Community (Similarity) Measures</font></u> (estimating various N-community similarity indices). The richness-based indices include the classic </br>
&nbsp;&nbsp;&nbsp;N-community Jaccard and Sørensen indices; the abundance-based indices include the Horn, Morisita-Horn, and the N-community Bray-Curtis indices.

<li> <u><font color="4464F1">    Part VI: Genetics (Differentiation) Measures</font></u>  (applying Part V to estimate allelic dissimilarity/differentiation among sub-populations based on</br> &nbsp;&nbsp;&nbsp;multiple-population genetics data).

NOTE: Part III of the original SPADE (Prediction) has been expanded to become an independent program known as iNEXT (iNterpolation and EXTrapolation) which provides diversity estimates for rarefied and extrapolated samples up to a maximum sample size or sample completeness specified by the user. The program iNEXT also features seamless plots of sample-size- and coverage-based rarefaction and extrapolation sampling curves. Users can download iNEXT from the same links given above.

The running procedures are summarized as follows. 


<p>   Step 1. Select an analysis part from the top menu of SpadeR window.<br>
<p>   Step 2. Select your data input format from the Data Setting on the left panel.<br>
           &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
           <font color="4464F1">(For Parts I and II, five data input formats are supported; for Parts III, IV and V, three data input formats are supported; for Part VI, only one format is supported.)</font><br>
           
<p>   Step 3. Check the &lt;Demo data&gt; radio button to load demo data (you can load your own data
           by checking &lt;Upload data&gt;).<br>
           &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
           <font color="4464F1">(To load your own incidence-raw data, you must specify the number of sampling units in each community in the left panel of SpadeR window.)</font><br>
           
<p>   Step 4. For Part V and VI, choose an order q (0, 1 or 2) and  the comparison  target (relative
           abundances or absolute abundances) you would<br>
           &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
           like to compute the (dis)similarity
           for any pair of samples.<br>
           
<p>   Step 5. Press the button &lt;Run!&gt; to get the output.<br>

<font color="4464F1">(We use either an analytic method or a bootstrap resampling method to compute s.e. and
confidence interval of an estimator. For the latter, the default number of bootstrap
replications is 100. You may specify a larger number to obtain more accurate results, but it
will take a longer time to get the output. Also, the bootstrap resampling procedures vary
with trial, meaning that two different runs for the same data may result in different s.e.
estimates and different confidence intervals.)</font>

Along the second row (output) menu, there are four output selection tabs:



<li> In the  **&quot;Estimation&quot;** tab panel, basic data information and various estimates are shown for the demo uploaded data.  You can click "download as</br>
&nbsp;&nbsp;&nbsp;txt file" at the bottom of the output to download these estimates as a file.</li>

<li> The functionality of the  **&quot;Visualization&quot;**  tab panel varies by part. For Species
and SharedSpecies parts, this tab panel shows various estimates </br>
&nbsp;&nbsp;&nbsp;and their confidence intervals to facilitate clear comparisons, while for Diversity Profile Estimation part, this tab panel shows the plots of</br> 
&nbsp;&nbsp;&nbsp;empirical diversity profile and estimated profile. These figures can be
downloaded by clicking "download as PNG file" at the bottom of the displayed<br>
&nbsp;&nbsp;&nbsp;figure.</li>

<li> The function of the **&quot;Data Viewer&quot;** tab panel is to display the first ten records of thedemo/uploaded data. The entire demo/uploaded data can be</br> 
&nbsp;&nbsp;&nbsp;downloaded by clicking "download as txt file" at the bottom of the displayed data.</li>

<li> In the **&quot;Introduction&quot;** panel, users can view a brief introduction to SpadeR and a summary of the running procedures for each part of analysis.</li>

<li> In the **&quot;User Guide&quot;** panel, a link will direct users to this user guide.</li>

<br>
To gain familiarity with the program, we suggest that users first run the demo data sets
included in SpadeR and check the output with that given in the SpadeR User's Guide (link
given above). Part of the output for each example is also interpreted in the guide to help users
understand the numerical results. The formulas for estimators featured in SpadeR with relevant
references are provided in the Appendix of the User's Guide.

<font color="red"><i>Please do not use SpadeR in any commercial form or distribute it to other people. Instead,
potential users should access the program directly through the SpadeR webpage (see above). If
you publish your work based on results from SpadeR, please make references to the relevant
papers mentioned in each section of SpadeR User's Guide and also use the following reference to cite SpadeR: </i></font>


Chao, A., Ma, K. H., Hsieh, T. C. and Chiu, C. H. (2015) Online Program SpadeR
(<u>S</u>pecies-richness <u>P</u>rediction <u>A</u>nd <u>D</u>iversity <u>E</u>stimation in <u>R</u>). Program and User's Guide
published at http://chao.stat.nthu.edu.tw/wordpress/software_download/.

<h3 style = "font-weight: bold">--------------------------------------------------------------------------------------------------------------------------------------</h3>

We recommend the following recent papers for pertinent background on biodiversity measures
and statistical analyses. These papers can be directly downloaded from Anne Chao's website.

<p>Chao, A.,and Chiu, C. H. (2012). Estimation of species richness and shared species richness. In N. Balakrishnan (ed). <i>Methods and Applications of Statistics in the Atmospheric and Earth Sciences</i>. p.76-111, Wiley, New York.<font color="4464F1">(Background on species richness and shared species richness estimation)</font>

<p>Chao, A., and Chiu, C. H. (2016). Nonparametric estimation and comparison of species richness.<i>Wiley Online Reference in the Life Science</i>. In: eLS. John Wiley & Sons, Ltd: Chichester. DOI: 10.1002/9780470015902.a0026329.<font color="4464F1"> (Background on comparing species richness across communities)</font>

<p>Chao, A., and Chiu, C. H. (2016). Bridging the variance and diversity decomposition approaches to beta diversity via similarity and differentiation measures. <i>Methods in Ecology and Evolution</i>. Early view at http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12551/abstract.
<font color="4464F1">(A unified theoretical framework on similarity/differentiation measures)</font>

<p>Chao, A., Chiu, C. H. and Jost, L. (2014). Unifying species diversity, phylogenetic diversity,
functional diversity, and related similarity/differentiation measures through Hill
numbers. <i>Annual Reviews of Ecology, Evolution, and Systematics</i> 45: 297-324.
<font color="4464F1">(A unified theoretical framework on diversity measures)</font>

Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H., Colwell, R. K. and Ellison, A. M.(2014). Rarefaction and extrapolation with Hill numbers: a framework for sampling and estimation in species diversity studies. <i>Ecological Monographs</i> 84:45-67.
<font color="4464F1">(Background on comparing diversity measures across communities)</font>

Chao, A. and Jost, L. (2015). Estimating diversity and entropy profiles via discovery rates of new species. <i>Methods in Ecology and Evolution</i>, 6, 873-882.<font color="4464F1">(A unified approach to estimating diversity in a community based on incomplete samples)</font>

Chao, A., Wang, Y. T. and Jost, L. (2013). Entropy and the species accumulation curve: a novel
entropy estimator via discovery rates of new species. <i>Methods in Ecology and Evolution</i>, 4,1091-1100.
<font color="4464F1">(A nearly optimal estimator of Shannon entropy/diversity based on incomplete samples)</font>

