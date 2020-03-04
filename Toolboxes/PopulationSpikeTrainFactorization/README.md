Population Spike Train Factorization Toolbox for Matlab
=======================================================

Toolbox for factorizing population spike trains into firing patterns and
activation coefficients.

If you use this toolbox, please cite the paper:
Onken, A., Liu, J. K., Karunasekara, P. P. C. R., Delis, I., Gollisch, T.,
Panzeri, S. (2016) Using matrix and tensor factorizations for the single-
trial analysis of population spike trains, PLoS Computational Biology.


Description
-----------

Advances in neuronal recording techniques are leading to ever larger
numbers of simultaneously monitored neurons. This poses the important
analytical challenge of how to capture compactly all sensory information
that neural population codes carry in their spatial dimension (differences
in stimulus tuning across neurons at different locations), in their
temporal dimension (temporal neural response variations), or in their
combination (temporally coordinated neural population firing). Here we
provide an implementation for tensor factorizations of population spike
trains along space and time. These factorizations decompose a dataset of
single-trial population spike trains into spatial firing patterns
(combinations of neurons firing together), temporal firing patterns
(temporal activation of these groups of neurons) and trial-dependent
activation coefficients (strength of recruitment of such neural patterns
on each trial). We also provide various additional functions for
obtaining the optimal number of patterns, for destroying information
in space or time, for extracted the first spikes of each trial and
each neuron and for calculating similarity of patterns.


Functions
---------

Run the script 'demo.m' to see an application of spatiotemporal NMF
and space-by-time NMF to the example data file 'example_data.mat'.
The toolbox contains the following functions:
* nmf              - Non-negative matrix factorization
* stnmf            - Spatiotemporal NMF
* select_n_m       - Select number of spatiotemporal modules
* sbtnmf           - Space-by-time NMF
* select_n_tm_n_sm - Select numbers of temporal and spatial modules
* ldacc            - LDA classification
* rank_order_train - Training function for the rank order decoder
* rank_order_test  - Prediction function for the rank order decoder
* count_spikes     - Transforms spike trains to spike count matrices
* latency_only     - Extract first spikes
* shuffle_space    - Shuffle population responses in space
* shuffle_time     - Shuffle population responses in time
* overlap          - Overlap between modules
* moddist          - Geodesic distance between two sets of modules


Related Methods
---------------

Other methods compared against in the paper are:
* Orthogonal Tucker-2 decomposition:
  Implemented in the tucker function of the N-way toolbox available
  at http://www.models.life.ku.dk/nwaytoolbox/
* Bayes Poisson Factor:
  Implemented in the Bayes Poisson Factor Toolbox for Matlab available
  at https://github.com/ch237/BayesPoissonFactor/
* Principal Component Analysis:
  Implemented in the princomp function of the Matlab Statistics and
  Machine Learning Toolbox
  https://uk.mathworks.com/help/stats/princomp.html
* Independent Component Analysis:
  Implemented in the FastICA Toolbox for Matlab available at
  https://research.ics.aalto.fi/ica/fastica/
* Factor Analysis:
  Implemented as Fast Factor Analysis for Matlab available at
  http://mlg.eng.cam.ac.uk/zoubin/software.html


License
-------

```text
The Population Spike Train Factorization Toolbox is free software; you
can redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation; either
version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, see <http://www.gnu.org/licenses/>.
```
