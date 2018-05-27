# TODO

- Integrate with ROS
    - input scans, output scans, processing progress
- Reduce looping while matching scans
- Handle edge cases in scan matching
    - occluded / exposed regions
- Scoring should take into account errors from distances - use probabilities
- Test with more sophisticated world examples
- Simulate robot going through path with error


# Done
- Performance:
  remove atan2 from difference: Keep scans in x,y format, use x/y or y/x for sorting angles 
- Vectorization for performance
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
