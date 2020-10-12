# ModelReconstruction
> Analysis of time series from stochastic processes.

# What's still to do

- [ ] Correctness of transform function
    
    - [x] ✅ ```RESOLVED``` It may be incorrect: s. artifacts and index errors. But these seem to decrease with increasing bins.
    > I'll ignore that for now.

    - [ ] 🟥 ```TODO``` On Lorenz attractor it works perfectly fine, when it comes to the non-noised. when it's noised the function is misbehaving.

- [x] ✅ ```RESOLVED``` n-dimensional implementation
    > It's gorgeous!

- [x] ✅ ```RESOLVED``` Van der Pol example

- [ ] 🟥 ```TODO``` Implement abstract function to analyze every kind of defined function ```f()```, so that the notebook looses weight, and not every kind of funtion needs a seperate section.

- [x] ✅ ```RESOLVED``` The calculated drift and diffusion coefficients need to be multiplied by some factor to match the analytical ones. This is b/c in the calculation we derivate by 

    ![formula](https://render.githubusercontent.com/render/math?math={\frac{1}{\tau}\text{%20where%20}\tau=1})

    But this needs to be projected to the "actual with" of that segment, aka ```np.linspace```.

- [x] ✅ ```RESOLVED``` Implement D^(2) function. How to calculate D^(2) with n-dimensions? How many D^(2)'s are in n-dim?

    > Done.

- [ ] 🟥 ```TODO``` Create class ```StochasticAnalysis```:
    
    ```.D_1()``` to get drift

    ```.D_2()``` to get diffusion

    ```.analyze()``` to analyze abstract and return graphs and stuff

- [x] ✅ ```RESOLVED``` Function to create time series with
    > I used ```scipy.integrate.solve_ivp()``` together with custom ```f(t, x)```.

- [x] ✅ ```RESOLVED``` How to solve FPE in python? Or how to implement it: A field and a density function?

    ![formula](https://render.githubusercontent.com/render/math?math={\frac{\partial}{\partial%20t}p(\vec{x},t%2B\tau|\vec{x},t)=(-\sum_i\frac{\partial}{\partial%20x_i}D_i^{(1)}(\vec{x},t)%2B\sum_{ij}%20+%20\frac{\partial}{\partial%20x_ix_j}D_ij^{(2)}(\vec{x},t))\cdot%20p(\vec{x},t%2B\tau|\vec{x},t)})

    > Euler integration. With IVP as one bin or more on phaseplot are one. the rest 0.

<!-- $$
    \frac{\partial}{\partial t} p(\vec{x}, t + \tau | \vec{x}, t) = (- \sum_i \frac{\partial}{\partial x_i} D_i^{(1)}(\vec{x}, t) + \sum_{ij} \frac{\partial}{\partial x_i x_j} D_ij^{(2)}(\vec{x}, t)) \cdot p(\vec{x}, t + \tau | \vec{x}, t)
$$ -->

- [ ] 🟥 ```TODO``` Solve FPE via euler integration.

- [x] ✅ ```RESOLVED``` How to actually derivate a field? 

    ![formula](https://render.githubusercontent.com/render/math?math={\frac{\partial}{\partial%20x_i}D_i^{(1)}\text{and}\frac{\partial}{\partial%20x_ix_j}D_ij^{(2)}\text{[1](2)}})

    > ```np.diff()```

<!-- $$
    \frac{\partial}{\partial x_i} D_i^{(1)} \text{ and } \frac{\partial}{\partial x_i x_j} D_ij^{(2)} \text{ [1](2)}
$$ -->

- [x] ✅ ```RESOLVED``` Is it true that 

    ![formula](https://render.githubusercontent.com/render/math?math={<a%20\cdot%20b>=<a>%20\cdot<b>})

    If it's true, it'll make calculation of D^(2) easier. See example of HO in ipynb.

    > It must be wrong! D^(2) needs to be calculated otherwise!

- [ ] 🟥 ```TODO``` Documentation as MD in ipynb.

- [ ] 🟥 ```TODO``` There is a mistake in D^(2) calculation. Most probably k and j need to be swapped.
    ```python
    for j in range(dimension):
        for k in range(dimension):
            a_grid[j][c] += (series[j][i + tau] - series[j][i]) * (series[k][i + tau] - series[k][i])
    ```