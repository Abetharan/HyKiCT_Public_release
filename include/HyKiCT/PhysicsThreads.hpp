//   HyKiCT - 1D Lagrangian Radiation-Hydrodyanmics code for testing Coupling with Kinetic codes.
//   Copyright (C) 2020- Abetharan Antony
//
//   This program is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <https://www.gnu.org/licenses/>. 
#ifndef RADIATIONTHREADS_HPP
#define RADIATIONTHREADS_HPP
#include <mutex>
#include <vector>
#include <condition_variable>
#include <petscksp.h>
#include <boost/asio/io_service.hpp>
#include <boost/bind/bind.hpp>
#include <boost/thread/thread.hpp>
#include "MultiRadTrans.h"
#include "GridData.h"
#include "Switches.h"
#include "FixedData.h"
#include "MatrixSolver.h"
#include "TimeStep.h"
#include "SNBSource.h"

class PhysicsThreads
{
    private:
        std::mutex mutex;
        std::condition_variable condition;
        int workCount;
        boost::asio::io_service ioService;
        boost::thread_group threadpool;
        boost::asio::io_service::work *work;
    public:
        PhysicsThreads(int np)
        {
            for(int p = 0; p < np; p++)
            {
                threadpool.create_thread(
                    boost::bind(&boost::asio::io_service::run, &ioService)
                );
            }
        }
        PhysicsThreads(){
            work = new boost::asio::io_service::work(ioService);
        }
        ~PhysicsThreads()
        {
            delete work;
            ioService.stop();
            try
            {
            threadpool.join_all();
            }
            catch(const std::exception) {}
        }
        void resize(int np)
        {
            for(int p = 0; p < np; p++)
            {
                threadpool.create_thread(
                    boost::bind(&boost::asio::io_service::run, &ioService)
                );
            }
        }
        void radiationLoopTransport(MultiRadTrans &radiationTransport, 
                                    GridData const *gridDataT, FixedData const &fixedData, 
                                    TimeStep const& timeData, Switches const&switchData)
        {
            if(switchData.RadFFOn)
            {
                radiationTransport.updateRadiationQuantities(gridDataT, fixedData, timeData, 1);
            }            
            if(switchData.RadFicksLaw)
            {
                radiationTransport.updateRadiationQuantities(gridDataT, fixedData, timeData, 2);
            }
            std::unique_lock<std::mutex> lock(mutex);
            workCount--;
            condition.notify_one();
        }
        void snbLoopTransport(SNBSource &snbTransport, 
                                    GridData const *gridDataT, FixedData const &fixedData)
        {
            snbTransport.snbHeatFlowCorrection(gridDataT, fixedData);
            std::unique_lock<std::mutex> lock(mutex);
            workCount--;
            condition.notify_one();
        }
        int update(std::vector<MultiRadTrans> &radiationTransport, 
                    GridData const *gridDataT, FixedData const &fixedData, 
                    TimeStep const& timeData, Switches const&switchData)
        {
            int j = 0;
            workCount = fixedData.RadNg;
            for(auto &i:radiationTransport)
            {
            pushJobs(i,  gridDataT, 
                    fixedData, timeData, switchData);
                j++;
            }
            {
                std::unique_lock<std::mutex> lock(mutex);
                condition.wait(lock, [&]{return workCount ==0;});
            }
            return 0;
        }

        int update(std::vector<SNBSource> &snbTransport,
                    GridData const *gridDataT, FixedData const &fixedData)
        {
            int j = 0;
            workCount = fixedData.Ng;
            for(auto &i:snbTransport)
            {
            pushJobs(i, gridDataT, 
                    fixedData);
                j++;
            }
            {
                std::unique_lock<std::mutex> lock(mutex);
                condition.wait(lock, [&]{return workCount ==0;});
            }
            return 0;
        }
        void pushJobs(MultiRadTrans &radiationTransport, 
                                    GridData const *gridDataT, FixedData const &fixedData, 
                                    TimeStep const& timeData, Switches const&switchData)
        {
            ioService.post(boost::bind(&PhysicsThreads::radiationLoopTransport,this, 
                radiationTransport, gridDataT, 
                fixedData, timeData, switchData));
        }
        void pushJobs(SNBSource &snbTransport, 
                                    GridData const *gridDataT, FixedData const &fixedData)
        {
            ioService.post(boost::bind(&PhysicsThreads::snbLoopTransport,this, 
                snbTransport,  gridDataT, 
                fixedData));
        }
};
#endif