# NeuralPSP
Benchmark of online Principal Subspace Projection algorithms along with a set of algorithms efficiently implemented in MATLAB

For the more complete Python version please go to the link [online-psp](https://github.com/flatironinstitute/online-psp)

## Installation

Clone the repository or unzip the source and add recursively folders from the src folder to the MATLAB path

## EXAMPLES 

### Basic Example

``` Matlab
k -> subspace dimension
d -> number of features
% we suggest to standardize data using the standardize_data function
[X,~,~] = standardize_data(X,0,0);

fsm = FSM(k, d, [], [], [], []);
for i = 1:n    
    fsm.fit_next(x(:,i)');
end

components = fsm.get_components([]);
```



### Detailed Example
For more detailed examples look at the demo_XXX.m files


## References
[1] Pehlevan, Cengiz, Anirvan M. Sengupta, and Dmitri B. Chklovskii. "Why do similarity matching objectives lead to Hebbian/anti-Hebbian networks?." Neural computation 30, no. 1 (2018): 84-124.

[2] Cardot, Hervé, and David Degras. "Online Principal Component Analysis in High Dimension: Which Algorithm to Choose?." arXiv preprint arXiv:1511.03688 (2015).
