cc = g++
prom = run
deps = $(shell find ./ -name "*.h")
src = $(shell find ./ -name "*.cpp")
obj = $(src:%.cpp=%.o)


$(prom):$(obj)
    $(cc) -fopenmp -o $(prom) $(obj)


%.o:%.cpp $(deps)
    $(cc) -std=c++11 -fopenmp -c $< -o $@


clean:
    rm -rf $(obj) $(prom)