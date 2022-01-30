
#include <cstdlib>
#include <cstring>
#include "lab.h"
#include <iostream>
#include "zhalninrv.h"
#include "akaykinsv.h"
#include "akimovada.h"
#include "artamonovav.h"
#include "bugreevaam.h"
#include "venediktovayap.h"
#include "vecherskiymp.h"
#include "denisovrv.h"
#include "zinkinakv.h"
#include "kaderovro.h"
#include "kochetkovpa.h"
#include "makarovaayu.h"
#include "melkonyanma.h"
#include "melyakinev.h"
#include "negryame.h"
#include "nikishkinev.h"
#include "nuyanzinma.h"
#include "pomelovaas.h"
#include "prokopenkoas.h"
#include "prokopenkods.h"
#include "rodkinav.h"
#include "ryabikinks.h"
#include "timovkinayu.h"
#include "turaevdv.h"
#include "fedinda.h"



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
  else if (strcmp(argv[1], "akaykinsv") == 0) {
    l = new akaykinsv();
  }
  else if (strcmp(argv[1], "akimovada") == 0) {
    l = new akimovada();
  }
  else if (strcmp(argv[1], "artamonovav") == 0) {
    l = new artamonovav();
  }
  else if (strcmp(argv[1], "bugreevaam") == 0) {
    l = new bugreevaam();
  }
  else if (strcmp(argv[1], "venediktovayap") == 0) {
    l = new venediktovayap();
  }
  else if (strcmp(argv[1], "vecherskiymp") == 0) {
    l = new vecherskiymp();
  }
  else if (strcmp(argv[1], "denisovrv") == 0) {
    l = new denisovrv();
  }
  else if (strcmp(argv[1], "zinkinakv") == 0) {
    l = new zinkinakv();
  }
  else if (strcmp(argv[1], "kaderovro") == 0) {
    l = new kaderovro();
  }
  else if (strcmp(argv[1], "kochetkovpa") == 0) {
    l = new kochetkovpa();
  }
  else if (strcmp(argv[1], "makarovaayu") == 0) {
    l = new makarovaayu();
  }
  else if (strcmp(argv[1], "melkonyanma") == 0) {
    l = new melkonyanma();
  }
  else if (strcmp(argv[1], "melyakinev") == 0) {
    l = new melyakinev();
  }
  else if (strcmp(argv[1], "negryame") == 0) {
    l = new negryame();
  }
  else if (strcmp(argv[1], "nikishkinev") == 0) {
    l = new nikishkinev();
  }
  else if (strcmp(argv[1], "nuyanzinma") == 0) {
    l = new nuyanzinma();
  }
  else if (strcmp(argv[1], "pomelovaas") == 0) {
    l = new pomelovaas();
  }
  else if (strcmp(argv[1], "prokopenkoas") == 0) {
    l = new prokopenkoas();
  }
  else if (strcmp(argv[1], "prokopenkods") == 0) {
    l = new prokopenkods();
  }
  else if (strcmp(argv[1], "rodkinav") == 0) {
    l = new rodkinav();
  }
  else if (strcmp(argv[1], "ryabikinks") == 0) {
    l = new ryabikinks();
  }
  else if (strcmp(argv[1], "timovkinayu") == 0) {
    l = new timovkinayu();
  }
  else if (strcmp(argv[1], "turaevdv") == 0) {
    l = new turaevdv();
  }
  else if (strcmp(argv[1], "fedinda") == 0) {
    l = new fedinda();
  }

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
  std::cout << "Usage:\n\n  " << name << " <fio> <lab_number>\n";
}
    
