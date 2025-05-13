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

