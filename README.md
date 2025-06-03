# Solar Panel Optimization System using Meta-Heuristics

## Project Overview
MATLAB-based optimization framework for residential solar panel configurations using four meta-heuristic algorithms. Balances:
- Energy production
- Installation costs
- Roof space utilization
- Thermal efficiency

## Algorithm Comparison

| Algorithm        | Type       | Runtime | Fitness | Strength          |
|------------------|------------|---------|---------|-------------------|
| Genetic (GA)     | Population | 7.5 min | 18.40   | Highest accuracy  |
| Grey Wolf (GWO)  | Population | 12 min  | 18.70   | Fast convergence  |
| Multi-Verse (MVO)| Population | 4 min   | 18.74   | Best balance      |
| Simulated Annealing | Trajectory | 10 min | 18.74   | Reliable baseline |


## Mathematical Model
### Objective Function
-
  ```math
    \min \left( \text{Cost} + \frac{A_0 - \sum(n_i \times A_i)}{A_0} + \sum\frac{(T_r-5)\times M_T}{m} - 0.9\times\sum(\text{Power}_i\times n_i) \right)
        
### Key Constraints
- **Power loss** ≤ 15%
- **Wasted area** ≤ 10%
- **Max 5 panel types**
- User-defined budget limit

## Implementation

### Input Parameters
- Roof dimensions (length × width)
- Panel specifications:
  - Type (e.g., Panasonic, Axitec)
  - Wattage (200W-300W)
  - Cost ($150-$220 per panel)
  - Area coverage (1.4-1.8 m²)
  - Temperature coefficients

### Optimization Process
1. Initialize population/solution
2. Evaluate against constraints
3. Apply algorithm-specific operators:
   - GA: Crossover and mutation
   - GWO: Alpha/beta/delta hierarchy
   - MVO: White/black hole cosmology
   - SA: Geometric cooling
4. Check convergence
5. Output best solution

## System Architecture
- 
    ```mermaid
      graph TD
          A[Input Parameters] --> B[Algorithm Selection]
          B --> C[Optimization Engine]
          C --> D[Constraint Validation]
          D --> E[Solution Refinement]
          E --> F[Output Generation]

## Results Summary
- **Best performer**: MVO (4min runtime, 18.74 fitness)
- **Most accurate**: GA (18.40 fitness)
- All solutions met:
  - <10% area waste
  - <15% power loss
  - ≤5 panel types limit

## Applications
- Residential solar planning
- Renewable energy studies
- Sustainable architecture

## Future Work
- 3D roof modeling integration
- Real-time weather data
- Mobile app interface
- Hybrid algorithms

*Developed by: Ahmed Mady, Mohammed Nada, Omar Shehata, Amro Abdellatif, Mostafa Kashaf, Marwan Sallam, Catherine Elias*# Solar Panel Optimization System using Meta-Heuristics

## Project Overview
MATLAB-based optimization framework for residential solar panel configurations using four meta-heuristic algorithms. Balances:
- Energy production
- Installation costs
- Roof space utilization
- Thermal efficiency

## Algorithm Comparison

| Algorithm        | Type       | Runtime | Fitness | Strength          |
|------------------|------------|---------|---------|-------------------|
| Genetic (GA)     | Population | 7.5 min | 18.40   | Highest accuracy  |
| Grey Wolf (GWO)  | Population | 12 min  | 18.70   | Fast convergence  |
| Multi-Verse (MVO)| Population | 4 min   | 18.74   | Best balance      |
| Simulated Annealing | Trajectory | 10 min | 18.74   | Reliable baseline |

## Key Constraints
- **Power loss** ≤ 15%
- **Wasted area** ≤ 10%
- **Max 5 panel types**
- User-defined budget limit

## Implementation

### Input Parameters
- Roof dimensions (length × width)
- Panel specifications:
  - Type (e.g., Panasonic, Axitec)
  - Wattage (200W-300W)
  - Cost ($150-$220 per panel)
  - Area coverage (1.4-1.8 m²)
  - Temperature coefficients

### Optimization Process
1. Initialize population/solution
2. Evaluate against constraints
3. Apply algorithm-specific operators:
   - GA: Crossover and mutation
   - GWO: Alpha/beta/delta hierarchy
   - MVO: White/black hole cosmology
   - SA: Geometric cooling
4. Check convergence
5. Output best solution

## Results Summary
- **Best performer**: MVO (4min runtime, 18.74 fitness)
- **Most accurate**: GA (18.40 fitness)
- All solutions met:
  - <10% area waste
  - <15% power loss
  - ≤5 panel types limit

## Applications
- Residential solar planning
- Renewable energy studies
- Sustainable architecture

## Future Work
- 3D roof modeling integration
- Real-time weather data
- Mobile app interface
- Hybrid algorithms

*Developed by: Ahmed Mady, Mohammed Nada, Omar Shehata, Amro Abdellatif, Mostafa Kashaf, Marwan Sallam, Catherine Elias*
