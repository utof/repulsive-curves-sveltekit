# QUESTIONS (and out-there ideas)

- ~~project constraints kinda works but how does it work~~

- what if find the subvertices that cause the most amount of energy (the intersecting ones) and kinda disable them to allow nonplanar graphs to not be crazy? 

### thoughts:
- ~~wtf even L2 is getting stuck. something is up with the gradients.~~
- ~~Remove THE LINE SEARCH~~
- [ ] wanna have a nice archive of the struckthrough stuff to see my progress and see how my understanding of the project changes





# todo

### for my own sake
- [ ] WHAT IS LINE SEARCH??
- [ ] understand exactly how barycenter works
  - [ ] same for constraint-jacobian making
  - [ ] 

### plugin system/modularity/microserviceness

- [ ] ability to define a constraint



### logic

- [ ] fix analytical derivative by giving only cpp + js files to gpt prompt

console log only every 1 sec or at the end of operations

1. pass b and b0 from somewhere else
2. store? or something for avoiding prop drilling

constraint to never exceed total length of the graph?
- [ ] api

### controls and config

- [ ] alpha, beta slider.
      - [ ] some auto way to determine parameters
- [ ] step-size parameter with saving TO the code from browser.
- [ ] dont regenerate graph after resaving the code
      - [ ] show dont show kernel

[ ] move controls out of the way

#### config 
- [ ] hardgradcap e.g. 1000 magnitude cap or 100% of the average or something or calculate % based on complexity of vertices and edges

### ui

- [ ] nicer controls ui
- [ ] logs/errors in the ui not only console (any lib for that?)


### dev experience

[ ] move all the controls (start stop animation, step) in the controls.svelte, not just the sliders.
[ ] apply perturbation function should be in utils or something
[ ] projectConstraints, enforceBarycenter need to be in constraint file
[ ] include subvertices energies condition inside the differential, not as a sep funciton
[ ] logging within the app like last energy, gradient magnitude, arrows of general translation of barycenter or something etc etc

#### generalizations 
[ ] define general graph in graphutils to not prevent modifying all the graphs when new feature is created.

[ ] define general gradient and specify specifications in other functions (precodnitioned, regular) with using the general func or smth.		           


[ ] somehow high-levely define the check of z-coord to avoid ugly check like these all the time:
```
		const dz = (p2[2] !== undefined || p1[2] !== undefined) ? 
		           ((p2[2] || 0) - (p1[2] || 0)) : undefined;
```

[ ] is there a way to define generally subvertices, curvature of subedges etc.? Like an api/class that will have properties e.g. drawInfo Position, relativeness of motion, extensionness of graphs.