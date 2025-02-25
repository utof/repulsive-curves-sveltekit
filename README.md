thoughts:

wtf even L2 is getting stuck. something is up with the gradients.
but wtf does barycenter do. i think it doesnt work.
Remove THE LINE SEARCH

# QUESTIONS:

WHAT IS LINE SEARCH??
project constraints kinda works but how does it work
# todo:

### logic

[ ] fix analytical by giving only cpp + js files. upd 2/25/25 what??

console log only every 1 sec or at the end of operations

1. pass b and b0 from somewhere else
2. store? or something for avoiding prop drilling

constraint to never exceed total length of the graph?
[ ] api

### controls

- [ ] alpha, beta slider.
      -- [ ] some auto way to determine parameters
- [ ] step-size parameter with saving TO the code from browser.
- [ ] dont regenerate graph after resaving the code
      [ ] show dont show kernel

[ ] move controls out of the way

### ui

nicer controls ui


### dev experience

[ ] move all the controls (start stop animation, step) in the controls.svelte, not just the sliders.
[ ] apply perturbation function should be in utils or something
[ ] projectConstraints, enforceBarycenter need to be in constraint file
[ ] define general graph in graphutils to not prevent modifying all the graphs when new feature is created.

[ ] is there a way to define generally subvertices, curvature of subedges etc.? Like an api/class that will have properties
e.g. drawInfo Position, relativeness of motion, extensionness of graphs.