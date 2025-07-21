Compare the disparity curve overlap between the "best" scenario (all nodes, no ace), vs. ace nodes - point estimates, vs. ace nodes - distributions (sample).

# Simulations:

make some treats trees with stop.rule = list(living = 50) and quiet extinction rate (e.g. bd.params = list(speciaition = 1, extinction = 0.7)) use null.error = 100 and replicate = 20 to get 20 trees with always 50 living species.

from these "true" trees (i.e. that have fossil AND node values), create three subsets of trees:
 1. just living species (i.e. drop.fossils)
 2. the "true" tree (i.e. keep everything)
 3. the "laggerstate fossilisation" tree (i.e. keep 50%)
 4. the "normal fossilisation" tree (i.e. keep N% whereN is from litterature, e.g. from Tom Smith's paper)

# The comparisons

1. How good is ace at recovering the values? compare the differences in a "distance space" (i.e. for discrete comparison in a MORD space and for continuous, comparisons in a euclidean space)

Remove node values, and then run ace, then compare ace values to "true" values to end up with two comparisons
    1.a ace.point.estimate vs. true node value: test = abs(ace - true) (this gives you the deviation from truth) 
    1.b ace.sample vs. true node value (i.e. in which sd is the true value) test = how many sd away from the sampled mean is the true value.

2. How does that affect the morphospace?
    
    2.a Ordinate the shit out of it and measure euclidean distance between true and estimated
    2.b Sub test: do ace within the trait space: basically ordinate the data from #simulations and then run ace on the ordinated data (i.e. post ordination ace) and measure the euclidean distance between true and estimated

3. How does that affect disparity?

    - Compare the curves using a rank envelope for the disparity through time for results of 1.a.b and 2.a.b


    COuld weight fossils pre-extinction during ace so that survivorship-bias is reduced?
Or alternatively, switch to a bayesian method of ace whereby the pre-extinction traitspace is used as a prior for the post-extinction estimates.
Do fossils help when mass extinction occurs?

Another idea:
- what is a better method: imputation of existing taxa sampling, or ancestral state estimation. 
- the rationale of this is that, often with these studies, fossils are degraded already with lots of missing data. 
- Is it better - do you recover true macroevolutionary patterns - if you use imputation methods of existing data, rather than ancestral state estimation.

Something to think about: precision vs accuracy

Could we assume that using slow rates is equivalent to using more broad coding system, and fast rates is small changes? Or should I link that to the number of coding states.

In the first, we converted Mk2 and BiSSE probabilities for each node into 0, 1, or ambiguous with this rule: if the probability was greater than 0.7, then the node state was assigned as 1, if the probability was less than 0.3, the state was assigned as 0, probabilities between 0.3 and 0.7 were considered as ambiguous. We then only considered “outright” errors, i.e. if the estimated and true states did not match. We refer to this as the quantised score. The second approach directly employed the node probabilities for Mk and BiSSE, and error scores of 0 (estimated and true values matched), 0.5 (ambiguous estimate) and 1 (estimated and true values were different) for MP.

GLM? to assess impact of parameters on probability of errors.

How can i ensure no negative values in ace - put a prior on the BM distribution??s


Do all nodes rather than just 49.
Clade extinct - simulates selective correlated extinction due to traits 

tom ezard for forams paper

So latest plots:
- disparity accuracy decreases as fossil sampling increases -> due to higher amount of nodes being estimated -> increases uncertainty. Wildly overestimates disparity, as uncertainty increases at these deeper nodes.

Not necessarily keeps uncertainty, but shows that when no state has likelihood above threshold, there is a lot of nodes that are left as 0/1, which makes differences between nodes more pronounced -> higher disparity.


Treats - do by max taxa. EXTINCTION hits at certain number of taxa

icons on poster.

old disparity figure for poster

Colours need to match
Use icons e.g. fossil sampling icon on x axis

For threshold method plot, show one scenario at 0.51, one at 0.9

For ace plot, do it with real animals -> wings and hands?


Change logo of relative threshold_method - needs just arrow pointing down to next in line - so that makes it 
 
PLot all 4 of the fossil sampling levels, with 1 headed arrow going back to show degradation - be explicit and show removal


BITS TO DO:
- make top panel slightly narrower to give more space to ace
- make it clear ace leads to disparity issues
- fast slow med needs to be highlighted in plots - DONE
- rewrite ace bit, make it clear that relative majority rule is more conservative, strict majority rule is permissive, write 100x matrices for sample
- perhaps takeaway the scaled likelihood bit - DONE
- write hamming as ancestral state error
- make it clear that y axis is good vs bad
- perhaps draw a line across 0
- make it clear that no majority rule propagates uncertainty.

on the bar charts, the state 0/1 needs to be larger.




- The actual ace plots - do they show accuracy change?? The `accuracy' of ace is the same across all 3, it just depends on how the likelihoods are interpreted. So I should perhaps look at doing a different way of interpreting that plot - how uncertainty increases with fossil sampling level -> this can be interpreted as average likelihood 

So I think follow Keating's error and uncertainty ideas -> see how error is affected and how uncertainty is propelled down pipeline.

SO even tho disparity estimates may appear wrong - they recover true pattern across time.


Done uncertainty - so with fossil sampling, uncertainty does not necessarily decrease (see high rates) because increased data increases uncertainty but does not necessarily lead to higher accuracy.
Can see this with the strict majority rule - the increase in fossil sampling increases accuracy but decreases uncertainty - lower uncertainty does not necessarily mean higher accuracy.


Should I change these plots to actual matrix preservation -> i.e fossil degredation, might be a more interesting parameter than fossil sampling.
Because obviously they will use the max fossils they can, but i think it would be useful to see where ancestral state estimation becomes more useful than just using tip data at what level of matrix preservation.



