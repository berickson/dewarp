# TODO
- Integrate with ROS
    - input scans, output scans, processing progress
- Investigate vectorization for performance
- Better hill climbing
    - change twiddle into hill-climb
    - try simulated annealing
    - idea - why not have a "leave this hill" in hill climbing that goes to other side of minimum? Maybe it only works in 1d? 


# Done
- Basic scan matching
- untwisting distortions from robot motion
- Make scan matching faster
    - Fast algorithm to move scan to a proposed pose
    - Each point can be changed to (r,theta) in new pose, 
    - compare points of two scans in order of theta
        - Use interpolation
        - tbd: what to do if theta gets out of order? (occlusion)
            - probably throw away further points by nulling
        - tbd: what to do about inconsistencies in points of view, like uncovered scenary? ()
