# Online-PSP
Efficient MATLAB implementation of online Principal Subspace Projection algorithms (Fast Similarity Matching[1], Incremental PCA[2,3], and Candid Covariance Incremental PCA[2,4])

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
For more detailed examples explore the demo_XXX.m files


## References
[1] Pehlevan, Cengiz, Anirvan M. Sengupta, and Dmitri B. Chklovskii. "Why do similarity matching objectives lead to Hebbian/anti-Hebbian networks?." Neural computation 30, no. 1 (2018): 84-124.

[2] Cardot, Hervé, and David Degras. "Online Principal Component Analysis in High Dimension: Which Algorithm to Choose?." arXiv preprint arXiv:1511.03688 (2015).

[3] Arora, R., Cotter, A., Livescu, K. and Srebro, N., 2012, October. Stochastic optimization for PCA and PLS. In Communication, Control, and Computing (Allerton), 2012 50th Annual Allerton Conference on (pp. 861-868). IEEE.

[4] Weng, J., Zhang, Y. and Hwang, W.S., 2003. Candid covariance-free incremental principal component analysis. IEEE Transactions on Pattern Analysis and Machine Intelligence, 25(8), pp.1034-1040.


## License
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
