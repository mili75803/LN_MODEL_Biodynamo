// -----------------------------------------------------------------------------
//
// Copyright (C) 2021 CERN & University of Surrey for the benefit of the
// BioDynaMo collaboration. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------
//
// A simulation of a conceptual model of a cancer tumor growth
//

#ifndef DEMO_TUMOR_CONCEPT_H_
#define DEMO_TUMOR_CONCEPT_H_

#include "biodynamo.h"
#include "core/util/random.h"
#include <chrono>
#include <mutex>

namespace bdm {
namespace tumor_concept {


class DC : public Cell {  // our object extends the Cell object
                              // create the header with our new data member
  BDM_AGENT_HEADER(DC, Cell, 1);

 public:
  DC() {}
  explicit DC(const Real3& position) : Base(position) {}
  virtual ~DC() {}

  // getter and setter for our new data member
  void SetCellColor(int cell_color) { cell_color_ = cell_color; }
  int GetCellColor() const { return cell_color_; }

 private:
  // declare new data member and define their type
  // private data can only be accessed by public function and not directly
  int cell_color_;

};


// Define my custom cell MyCell, which extends Cell by adding extra data
// members: cell_color and can_divide
class MyCell : public Cell {  // our object extends the Cell object
                              // create the header with our new data member
  BDM_AGENT_HEADER(MyCell, Cell, 1);

 public:
  MyCell() {}
  explicit MyCell(const Real3& position) : Base(position) {}
  virtual ~MyCell() {}


  void InitializeRandomGenerators(const uint64_t seedTheta, const uint64_t seedPhi) {
    randomTheta_ = new Random();
    randomTheta_->SetSeed(seedTheta);
    randomPhi_ = new Random();
    randomPhi_->SetSeed(seedPhi);
  }

  Random* GetRandomTheta() {return randomTheta_; }
  Random* GetRandomPhi() {return randomPhi_; }


  void SetCellColor(int cell_color) { cell_color_ = cell_color; }
  int GetCellColor() const { return cell_color_; }

  real_t max_bound, min_bound;
  int isDone;
  int time_step;
  std::string filename;
  int ID;
  real_t ntcell;
  real_t dis_x, dis_y, dis_z;
  
  std::vector<Real3>DC_Pos;
  void storeDCPos(std::vector<Real3>& DC_Pos_){DC_Pos = DC_Pos_;}
  std::vector<Real3>& getDCPos(){return DC_Pos;}
  
 private:
  // declare new data member and define their type
  // private data can only be accessed by public function and not directly
  int cell_color_;
  Random* randomTheta_;
  Random* randomPhi_;
};

// Define growth behaviour
struct Growth : public Behavior {
  BDM_BEHAVIOR_HEADER(Growth, Behavior, 1);

  Growth() { AlwaysCopyToNew(); }
  virtual ~Growth() {}

  void Run(Agent* agent) override {
    real_t scale = 0.097*100;
    if (auto* cell = dynamic_cast<MyCell*>(agent)) {
      Real3 curr_position = cell->GetPosition();
      while(true ){ 
        real_t x_pos = curr_position[0];
        real_t y_pos = curr_position[1];
        real_t z_pos = curr_position[2];

        real_t th = cell->GetRandomTheta()->Uniform(0, 1);
        real_t theta = 2 * M_PI * th;
        real_t ph = cell->GetRandomPhi()->Uniform(0, 1);
        
        real_t phi = acos(1 - 2 * ph);
        real_t new_step_x = scale * sin(phi) * cos(theta);
        real_t new_step_y = scale * sin(phi) * sin(theta);
        real_t new_step_z = scale * cos(phi);

        x_pos = x_pos + new_step_x;
        y_pos = y_pos + new_step_y;
        z_pos = z_pos + new_step_z;

      
        if((x_pos >= cell->min_bound && x_pos <= cell->max_bound) && (y_pos >= cell->min_bound && y_pos <= cell->max_bound) && (z_pos >= cell->min_bound && z_pos <= cell->max_bound)){
          Real3 cell_movements = {x_pos, y_pos, z_pos};
             // cell->positionTC_ = cell_movements;
          cell->SetPosition(cell_movements);  
          cell->time_step++;
          break;
        }
      }
    }
  }
};

inline int Simulate(int argc, const char** argv) {
  auto set_param = [](Param* param) {
    param->bound_space = Param::BoundSpaceMode::kClosed;
    param->min_bound = -1000;
    param->max_bound = 1000;  // cube of 100*100*100
  };

  Simulation simulation(argc, argv, set_param);
  auto* rm = simulation.GetResourceManager();
  auto* param = simulation.GetParam();
  auto* myrand = simulation.GetRandom();

  size_t nb_of_cells = 100000;  // number of cells in the simulation
  real_t x_coord, y_coord, z_coord;

  for (size_t i = 0; i < nb_of_cells; ++i) {
    // random real_t between maximumm boundary and minimum boundary of the cube (LN)
    x_coord = 0;
    y_coord = 500;
    z_coord = 500;

    uint64_t seedTheta = std::chrono::system_clock::now().time_since_epoch().count();
    uint64_t seedPhi = std::chrono::system_clock::now().time_since_epoch().count();

    // creating the cell at position x, y, z and setting the cell parameters
    MyCell* cell = new MyCell({x_coord, y_coord, z_coord});
    cell->InitializeRandomGenerators(seedTheta, seedPhi);
    cell->SetDiameter(10);
    cell->ntcell = nb_of_cells;
    cell->isDone = 0;
    cell->max_bound = param->max_bound;
    cell->min_bound = param->min_bound;
    cell->time_step = 0;
    cell->ID = i+1;
    cell->AddBehavior(new Growth());
    // will vary from 0 to 5. so 6 different layers depending on y_coord
    cell->SetCellColor(static_cast<int>((y_coord / param->max_bound * 6)));
    rm->AddAgent(cell);  // put the created cell in our cells structure
  }

  // Run simulation
  simulation.GetScheduler()->Simulate(500);

  std::cout << "Simulation completed successfully!" << std::endl;
  return 0;
}

}  // namespace tumor_concept
}  // namespace bdm

#endif  // DEMO_TUMOR_CONCEPT_H_
