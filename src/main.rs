use std::env;
use rand::Rng;
use std::time::Instant;
use std::cmp::Ordering;
use std::fmt::Formatter;
use std::fmt::Display;
use std::fmt::Error;

//functions to generate random numbers
fn percentage() -> f32 {
    rand::random()
}

fn intiger(min: usize, max: usize) -> usize {
    rand::thread_rng().gen_range(min..max)
}

fn float64(min: f64, max: f64) -> f64 {
    rand::thread_rng().gen_range(min..max)
}

//data struct for the entire population
pub struct AgentPool {
    agents: Vec<Agent>,
    mutation_prob: f32,
    fitness_value: f32,
    generation: usize,
    max_gen: usize,
    num_elite: usize,
    crossover_probability: f32,
}

//implimentation of the pool
impl AgentPool {
    //function to create the first generation
    pub fn initialize(num_agents: usize, mutation_prob: f32, genome_size: usize, min_value: f64, max_value: f64) -> Self {
        let agents = (0..num_agents).map(|_| Agent::random_genome(genome_size, min_value, max_value)).collect();
        Self {
            fitness_value: 1.0,// system value, starts at 1 to represent maximum error
            generation: 0,// system value, itteration counter
            agents,
            mutation_prob,
            max_gen: 10000,// system value, maximum program itteration limit
            num_elite: 5,// number of elites in the population (not a ratio)
            crossover_probability: 0.01 // probability of chiasmata x2 for left right
        }
    }
    //function to rank the populations fitness useing the index as a key. index 0 is the fittest
    pub fn order(&mut self,) {
        self.agents.iter_mut().for_each(|a| a.fitness_function());
        self.agents.sort_by(|a, b| {
            let o_a = a.fitness();
            let o_b = b.fitness();
            o_a.partial_cmp(&o_b).unwrap_or(Ordering::Equal)
        });
        self.fitness_value = self.agents[0].fitness();
    }
    //function to create the next generation
    pub fn next_gen(&mut self, min_value: f64, max_value: f64) {
        let n = self.agents.len();
        let half_len = n / 2;
        for i in 0..half_len + n % 2 {
            if i < self.num_elite{ //section for elite processing
                let a = &self.agents[i];
                let b = &self.agents[intiger(0, half_len)];
                let mut solo = a.fuse(b);

                if percentage() < self.crossover_probability {
                    solo.crossover_left(b);
                }
                if percentage() < self.crossover_probability {
                    solo.crossover_right(b);
                }
                if percentage() < self.mutation_prob {
                    solo.point_mutation(min_value, max_value);
                }
            self.agents[intiger(self.num_elite, n)] = solo;
            }
            else // section for regular processing
            {
                let a = &self.agents[i];
                let b = &self.agents[i + 1];
                let mut solo = a.fuse(b);

                if percentage() < self.crossover_probability {
                    solo.crossover_left(b);
                }
                if percentage() < self.crossover_probability {
                    solo.crossover_right(b);
                }
                if percentage() < self.mutation_prob {
                    solo.point_mutation(min_value, max_value);
                }
                self.agents[half_len + i] = solo;
            }
    }
    self.generation += 1;
    }
    //functions to access values
    pub fn best_fitness(&self) -> f32 {
        self.fitness_value}

    pub fn generation(&self) -> usize {
        self.generation}

    pub fn fittest(&self) -> Option<&Agent> {
        self.agents.get(0)}
}

//data struct for individual agents
pub struct Agent {
    genes: Vec<f64>,
    fitness: f32
}

// agent implimentation
impl Agent {
    // function to create a fully randomised genome
    pub fn random_genome(genome_size: usize, min_value: f64, max_value: f64) -> Self {
        let genes = (0..genome_size).map(|_| float64(min_value, max_value)).collect();
        Self {
            fitness: 0.0,
            genes,
        }
    }
    //function to fuse two genomes with random compoments from two parent genomes
    pub fn fuse(&self, other: &Self) -> Self {
        let num_genes = self.genes.len();
        let genes = (0..num_genes).map(|i| {
                if intiger(0, 1) == 1 {
                    self.genes[i]
                } else {
                    other.genes[i]}}
                ).collect();

        Self {
            genes,
            fitness: 0.0,
        }
    }
    //crossiver function to overwrite one section of a genome from left to center with another genomes vlaues
    pub fn crossover_left(&self, other: &Self) -> Self {
        let num_genes = self.genes.len();
        let crossover_point = intiger(0, num_genes/2);
        let genes = (0..num_genes).map(|i| {
                if i < crossover_point {
                    self.genes[i]
                } else {
                    other.genes[i]
                }
            }).collect();

        Self {
            genes,
            fitness: 0.0,
        }
    }
    //crossiver function to overwrite one section of a genome from right to center with another genomes vlaues
    pub fn crossover_right(&self, other: &Self) -> Self {
        let num_genes = self.genes.len();
        let crossover_point = intiger(0, num_genes/2);
        let genes = (0..num_genes).map(|i| {
                if i > crossover_point {
                    self.genes[i]
                } else {
                    other.genes[i]
                }
            }).collect();

        Self {
            genes,
            fitness: 0.0,
        }
    }
    //function to introduce random genes
    pub fn point_mutation(&mut self, min_value: f64, max_value: f64) {
        let index = intiger(0, self.genes.len());
        self.genes[index] = float64(min_value, max_value);
    }
    //function to dfine fitness (this must be customised by the user to achive a specific task)
    pub fn fitness_function(&mut self) {
        self.fitness = {
            let sum: f64 = self.genes.iter().sum();
            let sum: f64 = ((sum).exp() - (-sum).exp())/((sum).exp() + (-sum).exp());
            sum as f32
        }
    }

    pub fn fitness(&self) -> f32 {
        self.fitness
    }
}

//function to run the simulation and tie everythign together
pub fn simulation(num_agents: usize, mutation_prob: f32, genome_size: usize, min_value: f64, max_value: f64) {

    println!("Spawning Agents");
    let start_time = Instant::now();
    let mut s = AgentPool::initialize(num_agents, mutation_prob, genome_size, min_value, max_value);
    while s.best_fitness() > -0.1 {
        s.order();
        s.next_gen(min_value, max_value);
        if s.generation() == s.max_gen {
            println!("\n[!] ALERT: no satisfactory solution found\n(Maximum iteration limit exceeded)");
            break;
        }
        println!("generation: {}\t best fitness {}", s.generation(), s.best_fitness())
    }
    //conclusion report
    println!("\nEnviroment Terminated\n");
    println!("Generations:\t{}\n", s.generation());
    println!("Optimal Solution:\nscore({})\n{}\n", s.best_fitness(), s.fittest().unwrap());
    println!("Elapsed Time:\t{:?}\n", Instant::now() - start_time);
    println!("Exit Code: 0");
    
}
//function to allow the viewing of data
impl Display for Agent {
    fn fmt(&self, f: &mut Formatter) -> Result<(), Error> {
        let v = &self.genes;
        if v.len() == 0 {
            return Ok(());
        }
        for num in &v[0..v.len() - 1] {
            if let Err(e) = write!(f, "{}, ", &num.to_string()) {
                return Err(e);
            }
        }
        write!(f, "{}", &v[v.len() - 1])
    }
}

//main function
fn main() {
    let args: Vec<String> = env::args().collect();



    let num_agents = args //number of agents in the popuation
        .get(1)
        .and_then(|arg| Some(arg.parse().unwrap()))
        .unwrap_or(100);

    let mutation_prob = args //chance of a random mutation
        .get(2)
        .and_then(|arg| Some(arg.parse().unwrap()))
        .unwrap_or(0.5);

    let genome_size = args //number of genes contained by each agent (must account for fitness function operation)
        .get(3)
        .and_then(|arg| Some(arg.parse().unwrap()))
        .unwrap_or(10);

    let min_value = args // minimum value a gene can take
        .get(4)
        .and_then(|arg| Some(arg.parse().unwrap()))
        .unwrap_or(0.0);

    let max_value = args // maximum value a gene can take
        .get(5)
        .and_then(|arg| Some(arg.parse().unwrap()))
        .unwrap_or(1.0);


    simulation(num_agents, mutation_prob, genome_size, min_value, max_value);
}