cmake_minimum_required(VERSION 3.25)
project(neuroscience C)

set(CMAKE_C_STANDARD 11)

add_executable(
        neuroscience
        main.c
        pkg/numerical_calculation/ordinary_differential_equation.h
        pkg/numerical_calculation/ordinary_differential_equation.c
        pkg/neuron_model/hodgkin_huxley.c
        pkg/neuron_model/hodgkin_huxley.h
)
