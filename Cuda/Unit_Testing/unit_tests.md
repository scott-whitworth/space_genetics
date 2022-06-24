<h1>Outline of unit tests completed</h1>

<h2>testing_rank unit tests </h2>

<h3> dominationTest()</h3>
- dominationTest() looked at both giveRank() and giveDistance() from optimization.cu.and
- giveRank was tested to make sure it was giving correct ranks, which depended on the function dominationCheck() and based on the tests this function worked 
as expected with no changes needed to it from the previous summer's version.
- giveRank() was updated then tested to factor in errorStatus to the non-dominated sort
- giveRank() was tested to verify the indexing of the many vectors involved, and it has been confirmed that the vectors front, newFront, dominates, and dominatedByCount are all indexed and used correctly. 
- the ranks given have been verified to be correct based on tests that used nans, posDiffs and speedDiffs that were equal, and large sets with random posDiffs and speedDiffs
- giveDistance() has functioned as expected with no changes needed to the version at the start of the summer (other than the shift to vectors)
- the tests ran on giveDistance() used nans, equal posDiffs and speedDiffs and large randomized sets. The known sets were hand calculated and then compared to with the computed values and they matched.
<h3>CAtest()</h3>
- CAtest() examined how the anneal was changing based on the values for the posDiff and speedDiff and how fast it would change based on different functions
- the test itself uses two random and small values for posDiff and speedDiff and puts these through a loop that will decrease them each run. Each time it runs it will calculate the cost and then the annealing
- an excel document was later made that better outlines this process and it was found that the best annealing function for an impact mission was linear and ^4 for rendezvous.
<h3>newDuplicateTest()</h3>
- this test was made to test a new way to detect duplicates, which was comparing there posDiffs and speedDiffs to a certain tolerance. (the tolerance decided is 1e-14)
- the goal of this test was to make sure the original was not marked as a duplicate ever, but all of its clones were to be marked as duplicates
- This test also made sure the duplicate and parent vectors were filled correctly and no duplicates were parents.
- to verify these objectives, a small set was used that sometimes had many duplicates, or just one, and the results of how these duplicates were marked was printed
- NOTE: this test is currently outdated as it uses the adult struct to mark duplicates but we now use errorStatus