import os
import shutil



def create_file(name):

    hpp = name + ".h"
    cpp = name + ".cpp"
    shutil.copy('zhalninrv.h', hpp)
    shutil.copy('zhalninrv.cpp', cpp)

    with open(hpp, "rt") as file:
        x = file.read()
    with open(hpp, "wt") as file:
        x = x.replace("zhalninrv",name)
        file.write(x)
    
    with open(cpp, "rt") as file:
        x = file.read()
    with open(cpp, "wt") as file:
        x = x.replace("zhalninrv",name)
        file.write(x)


files = (
    'akaykinsv',
    'akimovada',
    'artamonovav',
    'bugreevaam',
    'venediktovayap',
    'vecherskiymp',
    'denisovrv',
    'zinkinakv',
    'kaderovro',
    'kochetkovpa',
    'makarovaayu',
    'melkonyanma',
    'melyakinev',
    'negryame',
    'nikishkinev',
    'nuyanzinma',
    'pomelovaas',
    'prokopenkoas',
    'prokopenkods',
    'rodkinav',
    'ryabikinks',
    'timovkinayu',
    'turaevdv',
    'fedinda'
)


for file in files:
    create_file(file)


with open("Makefile", "wt") as file:
    x = 'CC=g++\n'\
        'CFLAGS=-c -Wall -g\n'\
        'LDFLAGS= -g\n'\
        'SOURCES= \\\n'
    for line in files:
        x = x + line + ".cpp\\\n"

    x = x + 'zhalninrv.cpp\\\n'\
        'lab.cpp\\\n'\
        'main.cpp\n'

    x = x +  """

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=vvm

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
\t$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
\t$(CC) $(CFLAGS) $< -o $@

clean:
\trm *.o vvm

        """

    file.write(x)

with open("CMakeLists.txt", "wt") as file:
    x = """
cmake_minimum_required(VERSION 3.0)
project(vvm)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(vvm
"""
    for line in files:
        x = x + line + ".cpp "    

    x = x + '\nzhalninrv.cpp lab.cpp main.cpp)'

    file.write(x)


with open("main.cpp", "wt") as file:
    x = """
#include <cstdlib>
#include <cstring>
#include "lab.h"
#include <iostream>
#include "zhalninrv.h"
"""
    for name in files:
        x = x + '#include "'+name+'.h"\n'

    x = x + """


void print_usage(char* name);


int main(int argc, char** argv)
{
  if (argc < 3) {
    print_usage(argv[0]);
    return 0;
  }

  lab *l = NULL;
  if (strcmp(argv[1], "zhalninrv") == 0) {
    l = new zhalninrv();
  }
"""
    for name in files:
        x = x + '  else if (strcmp(argv[1], "'+name+'") == 0) {\n'\
                '    l = new '+name+'();\n'\
                '  }\n'
    
    x = x + """
  else  {
    print_usage(argv[0]);
    return 0;
  }

  l->read_file();
  l->run(atoi(argv[2]));
  l->write_result();
  l->check_result();

  //delete l; // TODO:
  return 0;
}


void print_usage(char* name)
{
  std::cout << "Usage:\\n\\n  " << name << " <fio> <lab_number>\\n";
}
    
"""
    file.write(x)