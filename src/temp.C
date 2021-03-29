              // Time-stamp: <2002-01-18 21:51:36 piet>
             //================================================================
            //                                                                |
           //           /__----__                         ........            |
          //       .           \                     ....:        :.          |
         //       :                 _\|/_         ..:                         |
        //       :                   /|\         :                     _\|/_  |
       //  ___   ___                  _____                      ___    /|\   |
      //  /     |   \    /\ \     / |   |   \  / |        /\    |   \         |
     //  |   __ |___/   /  \ \   /  |   |    \/  |       /  \   |___/         |
    //   |    | |  \   /____\ \ /   |   |    /   |      /____\  |   \     \/  |
   //     \___| |   \ /      \ V    |   |   /    |____ /      \ |___/     |   |
  //                                                                      /   |
 //              :                       _/|     :..                    |/    |
//                :..               ____/           :....          ..         |
/*   o   //          :.    _\|/_     /                   :........:           |
 *  O  `//\                 /|\                                               |
 *  |     /\                                                                  |
 *=============================================================================
 *
 *  nbody_sh1.C:  an N-body integrator with a shared but variable time step
 *                (the same for all particles but changing in time), using
 *                the Hermite integration scheme.
 *                        
 *                ref.: Hut, P., Makino, J. & McMillan, S., 1995,
 *                      Astrophysical Journal Letters 443, L93-L96.
 *                
 *  note: in this first version, all functions are included in one file,
 *        without any use of a special library or header files.
 *_____________________________________________________________________________
 *
 *  usage: nbody_sh1 [-h (for help)] [-d step_size_control_parameter]
 *                   [-e diagnostics_interval] [-o output_interval]
 *                   [-t total_duration] [-i (start output at t = 0)]
 *                   [-x (extra debugging diagnostics)]
 * 
 *         "step_size_control_parameter" is a coefficient determining the
 *            the size of the shared but variable time step for all particles
 *
 *         "diagnostics_interval" is the time between output of diagnostics,
 *            in the form of kinetic, potential, and total energy; with the
 *            -x option, a dump of the internal particle data is made as well
 * 
 *         "output_interval" is the time between successive snapshot outputs
 *
 *         "total_duration" is the integration time, until the program stops
 *
 *         Input and output are written from the standard i/o streams.  Since
 *         all options have sensible defaults, the simplest way to run the code
 *         is by only specifying the i/o files for the N-body snapshots:
 *
 *            nbody_sh1 < data.in > data.out
 *
 *         The diagnostics information will then appear on the screen.
 *         To capture the diagnostics information in a file, capture the
 *         standard error stream as follows:
 *
 *            (nbody_sh1 < data.in > data.out) >& data.err
 *
 *  Note: if any of the times specified in the -e, -o, or -t options are not an
 *        an integer multiple of "step", output will occur slightly later than
 *        predicted, after a full time step has been taken.  And even if they
 *        are integer multiples, round-off error may induce one extra step.
 *_____________________________________________________________________________
 *
 *  External data format:
 *
 *     The program expects input of a single snapshot of an N-body system,
 *     in the following format: the number of particles in the snapshot n;
 *     the time t; mass mi, position ri and velocity vi for each particle i,
 *     with position and velocity given through their three Cartesian
 *     coordinates, divided over separate lines as follows:
 *
 *                      n
 *                      t
 *                      m1 r1_x r1_y r1_z v1_x v1_y v1_z
 *                      m2 r2_x r2_y r2_z v2_x v2_y v2_z
 *                      ...
 *                      mn rn_x rn_y rn_z vn_x vn_y vn_z
 *
 *     Output of each snapshot is written according to the same format.
 *
 *  Internal data format:
 *
 *     The data for an N-body system is stored internally as a 1-dimensional
 *     array for the masses, and 2-dimensional arrays for the positions,
 *     velocities, accelerations and jerks of all particles.
 *_____________________________________________________________________________
 *
 *    version 1:  Jan 2002   Piet Hut, Jun Makino
 *_____________________________________________________________________________
 */

#include  <iostream>
#include  <math.h>
#include  <cmath>                          // to include sqrt(), etc.
#include  <cstdlib>                        // for atoi() and atof()
#include  <unistd.h>                       // for getopt()
#include  <fstream>                        // for file output
#include  <string>
using namespace std;

typedef double  real;                      // "real" as a general name for the
                                           // standard floating-point data type

const int NDIM = 3;                        // number of spatial dimensions

void correct_step(real pos[][NDIM], real vel[][NDIM], 
                  const real acc[][NDIM], const real jerk[][NDIM],
                  const real old_pos[][NDIM], const real old_vel[][NDIM], 
                  const real old_acc[][NDIM], const real old_jerk[][NDIM],
                  int n, real dt);
void evolve(const real mass[], const real soft[], real pos[][NDIM], real vel[][NDIM],
            int n, real & t, real dt_param, real dt_dia, real dt_out,
            real dt_tot, bool init_out, bool x_flag, bool mw_flag, bool df_flag, 
            const real vcirc, const real ds, real tstep, real dtmin, real dtmax, int tsign);
void evolve_step(const real mass[], const real soft[], real pos[][NDIM], real vel[][NDIM],
                 real acc[][NDIM], real jerk[][NDIM], int n, real & t,
                 real dt, real & epot, real & coll_time, bool mw_flag, bool df_flag, 
                 const real vcirc, const real ds);
void get_acc_jerk_pot_coll(const real mass[], const real soft[], const real pos[][NDIM],
                           const real vel[][NDIM], real acc[][NDIM],
                           real jerk[][NDIM], int n, real & epot,
                           real & coll_time, bool mw_flag, bool df_flag, const real vcirc,
                           const real ds);
void get_snapshot(real mass[], real soft[], real pos[][NDIM], real vel[][NDIM], int n);
void predict_step(real pos[][NDIM], real vel[][NDIM], 
                  const real acc[][NDIM], const real jerk[][NDIM],
                  int n, real dt);
void put_snapshot(const real mass[], const real soft[], const real pos[][NDIM],
                  const real vel[][NDIM], int n, real t,
                  const real acc[][NDIM], const real jerk[][NDIM], const int nput);
bool read_options(int argc, char *argv[], real & dt_param, real & dt_dia,
                  real & dt_out, real & dt_tot, bool & i_flag, bool & x_flag,
                  bool & mw_flag, bool & df_flag, real & tstep, real & dtmin,
                  real & dtmax, int & tsign);

void write_diagnostics(const real mass[], const real soft[], const real pos[][NDIM],
                       const real vel[][NDIM], const real acc[][NDIM],
                       const real jerk[][NDIM], int n, real t, real epot,
                       int nsteps, real & einit, bool init_flag,
                       bool x_flag);

real jfit(real x);

/*-----------------------------------------------------------------------------
 *  main  --  reads option values, reads a snapshot, and launches the
 *            integrator
 *-----------------------------------------------------------------------------
 */

int main(int argc, char *argv[])
{
    real  dt_param = 0.03;     // control parameter to determine time step size, dt_param>0
    real  dt_dia = 1;          // time interval between diagnostics output,      dt_dia>0
    real  dt_out = 1;          // time interval between output of snapshots,     dt_out>0
    real  dt_tot = 10;         // duration of the integration,                   dt_tot>0
    bool  init_out = false;    // if true: snapshot output with start at t = 0
                               //          with an echo of the input snapshot
    bool  x_flag = false;      // if true: extra debugging diagnostics output
    bool  mw_flag = true;      // if true: use milky way potential
    bool  df_flag = false;     // if true: dynamical friction due to mw halo
    real  tstep = 0.;          // timestep to use (input by user, positive)
    real  dtmin = 0.;          // minimum timestep (input by user, positive)
    real  dtmax = 1e10;       // maximum timestep (input by user, positive)
    int   tsign = -1;          // backwards or forwards integration

    if (! read_options(argc, argv, dt_param, dt_dia, dt_out, dt_tot, init_out,
                       x_flag, mw_flag, df_flag, tstep, dtmin, dtmax, tsign))
        return 1;                // halt criterion detected by read_options()

    int n;                       // N, number of particles in the N-body system
    cin >> n;

    real t;                      // time
    cin >> t;

    real vcirc;            // vcirc
    cin >> vcirc;

    real ds;           // ds
    cin >> ds;

    real * mass = new real[n];                  // masses for all particles
    real * soft = new real[n];                  // softening param for all particles
    real (* pos)[NDIM] = new real[n][NDIM];     // positions for all particles
    real (* vel)[NDIM] = new real[n][NDIM];     // velocities for all particles

    get_snapshot(mass, soft, pos, vel, n);

    evolve(mass, soft, pos, vel, n, t, dt_param, dt_dia, dt_out, dt_tot, init_out,
           x_flag, mw_flag, df_flag, vcirc, ds, tstep, dtmin, dtmax, tsign);

    delete[] mass;
    delete[] soft;
    delete[] pos;
    delete[] vel;
}

/*-----------------------------------------------------------------------------
 *  read_options  --  reads the command line options, and implements them.
 *
 *  note: when the help option -h is invoked, the return value is set to false,
 *        to prevent further execution of the main program; similarly, if an
 *        unknown option is used, the return value is set to false.
 *-----------------------------------------------------------------------------
 */

bool read_options(int argc, char *argv[], real & dt_param, real & dt_dia,
                  real & dt_out, real & dt_tot, bool & i_flag, bool & x_flag, bool & mw_flag,
                  bool & df_flag, real & tstep, real & dtmin, real & dtmax, int & tsign)
{
    int c;
    while ((c = getopt(argc, argv, "hd:e:o:t:s:l:c:ixmyf")) != -1)
        switch(c){
            case 'h': cerr << "usage: " << argv[0]
                           << " [-h (for help)]"
                           << " [-d step_size_control_parameter]\n"
                           << "         [-e diagnostics_interval]"
                           << " [-o output_interval]\n"
                           << "         [-t total_duration]"
                           << " [-i (start output at t = 0)]\n"
                           << "         [-x (extra debugging diagnostics)]"
                           << " [-m (turn off MW potential)]\n"
                           << "         [-y (turn on dynamical friction)]"
                           << " [-s (timestep, fixed, positive)]"
                           << "         [-l (minimum timestep, positive)]\n"
                           << " [-c (maximum timestep, positive)]"
                           << "         [-f (forward time)]\n"

                           << endl;
                      return false;         // execution should stop after help
            case 'd': dt_param = atof(optarg);
                      break;
            case 'e': dt_dia = atof(optarg);
                      break;
            case 'i': i_flag = true;
                      break;
            case 'o': dt_out = atof(optarg);
                      break;
            case 't': dt_tot = atof(optarg);
                      break;
            case 'x': x_flag = true;
                      break;
	    case 'm': mw_flag = false;
	              break;
	    case 'y': df_flag = true;
	              break;
  	    case 's': tstep = atof(optarg);
	              break;
	    case 'l': dtmin = atof(optarg);
	              break;
	    case 'c': dtmax = atof(optarg);
	              break;
	    case 'f': tsign = 1;
	              break;

            case '?': cerr << "usage: " << argv[0]
                           << " [-h (for help)]"
                           << " [-d step_size_control_parameter]\n"
                           << "         [-e diagnostics_interval]"
                           << " [-o output_interval]\n"
                           << "         [-t total_duration]"
                           << " [-i (start output at t = 0)]\n"
                           << "         [-x (extra debugging diagnostics)]"
                           << " [-m (turn off MW potential)]"
                           << endl;
                      return false;        // execution should stop after error
            }

    return true;                         // ready to continue program execution
}

/*-----------------------------------------------------------------------------
 *  get_snapshot  --  reads a single snapshot from the input stream cin.
 *
 *  note: in this implementation, only the particle data are read in, and it
 *        is left to the main program to first read particle number and time
 *-----------------------------------------------------------------------------
 */

void get_snapshot(real mass[], real soft[], real pos[][NDIM], real vel[][NDIM], int n)
{
    for (int i = 0; i < n ; i++){
        cin >> mass[i];                       // mass of particle i
        cin >> soft[i];                       // softening param of particle i
        for (int k = 0; k < NDIM; k++)
            cin >> pos[i][k];                 // position of particle i
        for (int k = 0; k < NDIM; k++)
            cin >> vel[i][k];                 // velocity of particle i
    }
}
    
/*-----------------------------------------------------------------------------
 *  put_snapshot  --  writes a single snapshot on the output stream cout.
 *
 *  note: unlike get_snapshot(), put_snapshot handles particle number and time
 *-----------------------------------------------------------------------------
 */

void put_snapshot(const real mass[], const real soft[], const real pos[][NDIM],
                  const real vel[][NDIM], int n, real t,
                  const real acc[][NDIM], const real jerk[][NDIM], const int nput)
{
    cout.precision(16);                       // full double precision

    // Nbodies > 50000, outputting to a file
    if (n > 50000) {

      // Making the filename
      char *stem = "gal";
      char filename[50];
      char str[7];
      sprintf (str, "%d", nput);
      strcat (filename,stem);
      strcat (filename,str);
      strcat (filename,".out");
      //puts (filename);
    
      ofstream file ( filename, ios::out );

      //cout << filename << endl;
    
      if (file.is_open())
        {
    
          file << n << endl;                        // N, total particle number
          file << t << endl;                        // current time
          for (int i = 0; i < n ; i++){
            file << mass[i] << ' ';                      // mass of particle i
            file << soft[i];                      // softening param of particle i
            for (int k = 0; k < NDIM; k++)
              file << ' ' << pos[i][k];         // position of particle i
            for (int k = 0; k < NDIM; k++)
              file << ' ' << vel[i][k];         // velocity of particle i
            file << endl;
          }
          file.close();
        }


    // Nbodies < 50000, outputting to stdout
    } else {

      cout << n << endl;                        // N, total particle number
      cout << t << endl;                        // current time
      for (int i = 0; i < n ; i++){
          cout << mass[i] << ' ';                      // mass of particle i
          cout << soft[i];                      // softening param of particle i
          for (int k = 0; k < NDIM; k++)
              cout << ' ' << pos[i][k];         // position of particle i
          for (int k = 0; k < NDIM; k++)
               cout << ' ' << vel[i][k];         // velocity of particle i
          //for (int k = 0; k < NDIM; k++)
          //    cout << ' ' << acc[i][k];         // acc of particle i
          //for (int k = 0; k < NDIM; k++)
          //    cout << ' ' << jerk[i][k];         // jerk of particle i
          cout << endl;
      }

    }

}
    
/*-----------------------------------------------------------------------------
 *  write_diagnostics  --  writes diagnostics on the error stream cerr:
 *                         current time; number of integration steps so far;
 *                         kinetic, potential, and total energy; absolute and
 *                         relative energy errors since the start of the run.
 *                         If x_flag (x for eXtra data) is true, all internal
 *                         data are dumped for each particle (mass, position,
 *                         velocity, acceleration, and jerk).
 *
 *  note: the kinetic energy is calculated here, while the potential energy is
 *        calculated in the function get_acc_jerk_pot_coll().
 *-----------------------------------------------------------------------------
 */

void write_diagnostics(const real mass[], const real soft[], const real pos[][NDIM],
                       const real vel[][NDIM], const real acc[][NDIM],
                       const real jerk[][NDIM], int n, real t, real epot,
                       int nsteps, real & einit, bool init_flag,
                       bool x_flag)
{
    real ekin = 0;                       // kinetic energy of the n-body system
    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++)
            ekin += 0.5 * mass[i] * vel[i][k] * vel[i][k];

    real etot = ekin + epot;             // total energy of the n-body system

    if (init_flag)                       // at first pass, pass the initial
        einit = etot;                    // energy back to the calling function

   // cerr << "at time t = " << t << " , after " << nsteps
   //      << " steps :\n  E_kin = " << ekin
   //      << " , E_pot = " << epot
   //      << " , E_tot = " << etot << endl;
   // cerr << "                "
   //      << "absolute energy error: E_tot - E_init = "
   //      << etot - einit << endl;
   // cerr << "                "
   //      << "relative energy error: (E_tot - E_init) / E_init = "
   //      << (etot - einit) / einit << endl;

    cerr << t << ' ' << etot << ' '
         << ekin << ' ' << epot << ' ' << nsteps << endl;

    if (x_flag){
        cerr << "  for debugging purposes, here is the internal data "
             << "representation:\n";
        for (int i = 0; i < n ; i++){
            cerr << "    internal data for particle " << i+1 << " : " << endl;
            cerr << "      ";
            cerr << mass[i] << ' ';
            cerr << soft[i] << ' ';
            for (int k = 0; k < NDIM; k++)
                cerr << ' ' << pos[i][k];
            for (int k = 0; k < NDIM; k++)
                cerr << ' ' << vel[i][k];
            for (int k = 0; k < NDIM; k++)
                cerr << ' ' << acc[i][k];
            for (int k = 0; k < NDIM; k++)
                cerr << ' ' << jerk[i][k];
            cerr << endl;
        }
    }
}
    
/*-----------------------------------------------------------------------------
 *  evolve  --  integrates an N-body system, for a total duration dt_tot.
 *              Snapshots are sent to the standard output stream once every
 *              time interval dt_out.  Diagnostics are sent to the standard
 *              error stream once every time interval dt_dia.
 *
 *  note: the integration time step, shared by all particles at any given time,
 *        is variable.  Before each integration step we use coll_time (short
 *        for collision time, an estimate of the time scale for any significant
 *        change in configuration to happen), multiplying it by dt_param (the
 *        accuracy parameter governing the size of dt in units of coll_time),
 *        to obtain the new time step size.
 *
 *  Before moving any particles, we start with an initial diagnostics output
 *  and snapshot output if desired.  In order to write the diagnostics, we
 *  first have to calculate the potential energy, with get_acc_jerk_pot_coll().
 *  That function also calculates accelerations, jerks, and an estimate for the
 *  collision time scale, all of which are needed before we can enter the main
 *  integration loop below.
 *       In the main loop, we take as many integration time steps as needed to
 *  reach the next output time, do the output required, and continue taking
 *  integration steps and invoking output this way until the final time is
 *  reached, which triggers a `break' to jump out of the infinite loop set up
 *  with `while(true)'.
 *-----------------------------------------------------------------------------
 */

void evolve(const real mass[], const real soft[], real pos[][NDIM], real vel[][NDIM],
            int n, real & t, real dt_param, real dt_dia, real dt_out,
            real dt_tot, bool init_out, bool x_flag, bool mw_flag, bool df_flag, 
            const real vcirc, const real ds, real tstep, real dtmin, real dtmax, int tsign)
{
    cerr << "# Starting a Hermite integration for a " << n
         << "-body system,\n#   from time t = " << t 
         << " with time step control parameter dt_param = " << dt_param
         << "  until time " << t - dt_tot 
         << " ,\n#   with diagnostics output interval dt_dia = "
         << dt_dia << ",\n#   and snapshot output interval dt_out = "
         << dt_out << "." << endl;

    real (* acc)[NDIM] = new real[n][NDIM];          // accelerations and jerks
    real (* jerk)[NDIM] = new real[n][NDIM];         // for all particles
    real epot;                        // potential energy of the n-body system
    real coll_time;                   // collision (close encounter) time scale

    get_acc_jerk_pot_coll(mass, soft, pos, vel, acc, jerk, n, epot, coll_time, mw_flag, df_flag, vcirc, ds);

    int nsteps = 0;               // number of integration time steps completed
    int nput = 0;                 // number of outputs
    real einit;                   // initial total energy of the system

    write_diagnostics(mass, soft, pos, vel, acc, jerk, n, t, epot, nsteps, einit,
                      true, x_flag);
    if (init_out)                                    // flag for initial output
        put_snapshot(mass, soft, pos, vel, n, t, acc, jerk, 0);

    real t_dia = t + tsign * dt_dia;           // next time for diagnostics output, t<0, dt_...>0
    real t_out = t + tsign * dt_out;           // next time for snapshot output
    real t_end = t + tsign * dt_tot;           // final time, to finish the integration

    while (true){
        while (tsign*t < tsign*t_dia && tsign*t < tsign*t_out && tsign*t < tsign*t_end){   // t<0
	    real dt = tsign * dt_param * coll_time;          // dt<0, dt_param, coll_time >0

            // cerr << dt << " " << dt_param << " " << coll_time << endl;

            if (tstep != 0.)
                dt = tsign*tstep;              // using user supplied timestep
	    if (tsign*dt < dtmin)
                dt = tsign*dtmin;              // dt < dtmin, using dtmin
	    if (tsign*dt > dtmax)
                dt = tsign*dtmax;              // dt > dtmax, using dtmax

            // cerr << dt << " " << dt_param << " " << coll_time << endl;

	    //dt = -0.0001;
            evolve_step(mass, soft, pos, vel, acc, jerk, n, t, dt, epot, coll_time, mw_flag, df_flag, vcirc, ds);
            nsteps++;
        }
        if (tsign*t >= tsign*t_dia){             // tsign*t>0
            write_diagnostics(mass, soft, pos, vel, acc, jerk, n, t, epot, nsteps,
                              einit, false, x_flag);
            t_dia = t_dia + tsign * dt_dia;
        } 
        if (tsign*t >= tsign*t_out){              // tsign*t>0
            nput++;
            put_snapshot(mass, soft, pos, vel, n, t, acc, jerk, nput);
            t_out = t_out + tsign * dt_out;
        }
        if (tsign*t >= tsign*t_end)               // tsign*t>0
            break;
    }

    delete[] acc;
    delete[] jerk;
}

/*-----------------------------------------------------------------------------
 *  evolve_step  --  takes one integration step for an N-body system, using the
 *                   Hermite algorithm.
 *-----------------------------------------------------------------------------
 */

void evolve_step(const real mass[], const real soft[], real pos[][NDIM], real vel[][NDIM],
                 real acc[][NDIM], real jerk[][NDIM], int n, real & t,
                 real dt, real & epot, real & coll_time, bool mw_flag, bool df_flag,
                 const real vcirc, const real ds)
{
    real (* old_pos)[NDIM] = new real[n][NDIM];
    real (* old_vel)[NDIM] = new real[n][NDIM];
    real (* old_acc)[NDIM] = new real[n][NDIM];
    real (* old_jerk)[NDIM] = new real[n][NDIM];

    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++){
          old_pos[i][k] = pos[i][k];
          old_vel[i][k] = vel[i][k];
          old_acc[i][k] = acc[i][k];
          old_jerk[i][k] = jerk[i][k];
        }

    predict_step(pos, vel, acc, jerk, n, dt);
    //cerr << t << endl;
    get_acc_jerk_pot_coll(mass, soft, pos, vel, acc, jerk, n, epot, coll_time, mw_flag, df_flag, vcirc, ds);
    correct_step(pos, vel, acc, jerk, old_pos, old_vel, old_acc, old_jerk,
                 n, dt);
    t += dt;                // going backwards, dt<0, t<0

    delete[] old_pos;
    delete[] old_vel;
    delete[] old_acc;
    delete[] old_jerk;
}

/*-----------------------------------------------------------------------------
 *  predict_step  --  takes the first approximation of one Hermite integration
 *                    step, advancing the positions and velocities through a
 *                    Taylor series development up to the order of the jerks.
 *-----------------------------------------------------------------------------
 */

void predict_step(real pos[][NDIM], real vel[][NDIM], 
                  const real acc[][NDIM], const real jerk[][NDIM],
                  int n, real dt)
{
    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++){
            pos[i][k] += vel[i][k]*dt + acc[i][k]*dt*dt/2
                                      + jerk[i][k]*dt*dt*dt/6;
            vel[i][k] += acc[i][k]*dt + jerk[i][k]*dt*dt/2;
        }
}

/*-----------------------------------------------------------------------------
 *  correct_step  --  takes one iteration to improve the new values of position
 *                    and velocities, effectively by using a higher-order
 *                    Taylor series constructed from the terms up to jerk at
 *                    the beginning and the end of the time step.
 *-----------------------------------------------------------------------------
 */

void correct_step(real pos[][NDIM], real vel[][NDIM], 
                  const real acc[][NDIM], const real jerk[][NDIM],
                  const real old_pos[][NDIM], const real old_vel[][NDIM], 
                  const real old_acc[][NDIM], const real old_jerk[][NDIM],
                  int n, real dt)
{
    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++){
            vel[i][k] = old_vel[i][k] + (old_acc[i][k] + acc[i][k])*dt/2
                                      + (old_jerk[i][k] - jerk[i][k])*dt*dt/12;
            pos[i][k] = old_pos[i][k] + (old_vel[i][k] + vel[i][k])*dt/2
                                      + (old_acc[i][k] - acc[i][k])*dt*dt/12;
        }
}

/*-----------------------------------------------------------------------------
 *  get_acc_jerk_pot_coll  --  calculates accelerations and jerks, and as side
 *                             effects also calculates potential energy and
 *                             the time scale coll_time for significant changes
 *                             in local configurations to occur.
 *                                                  __                     __
 *                                                 |          -->  -->       |
 *               M                           M     |           r  . v        |
 *   -->          j    -->       -->          j    | -->        ji   ji -->  |
 *    a   ==  --------  r    ;    j   ==  -------- |  v   - 3 ---------  r   |
 *     ji     |-->  |3   ji        ji     |-->  |3 |   ji      |-->  |2   ji |
 *            | r   |                     | r   |  |           | r   |       |
 *            |  ji |                     |  ji |  |__         |  ji |     __|
 *                             
 *  note: it would be cleaner to calculate potential energy and collision time
 *        in a separate function.  However, the current function is by far the
 *        most time consuming part of the whole program, with a double loop
 *        over all particles that is executed every time step.  Splitting off
 *        some of the work to another function would significantly increase
 *        the total computer time (by an amount close to a factor two).
 *
 *  We determine the values of all four quantities of interest by walking
 *  through the system in a double {i,j} loop.  The first three, acceleration,
 *  jerk, and potential energy, are calculated by adding successive terms;
 *  the last, the estimate for the collision time, is found by determining the 
 *  minimum value over all particle pairs and over the two choices of collision
 *  time, position/velocity and sqrt(position/acceleration), where position and
 *  velocity indicate their relative values between the two particles, while
 *  acceleration indicates their pairwise acceleration.  At the start, the
 *  first three quantities are set to zero, to prepare for accumulation, while
 *  the last one is set to a very large number, to prepare for minimization.
 *       The integration loops only over half of the pairs, with j > i, since
 *  the contributions to the acceleration and jerk of particle j on particle i
 *  is the same as those of particle i on particle j, apart from a minus sign
 *  and a different mass factor.
 *-----------------------------------------------------------------------------
 */

void get_acc_jerk_pot_coll(const real mass[], const real soft[], const real pos[][NDIM],
                           const real vel[][NDIM], real acc[][NDIM],
                           real jerk[][NDIM], int n, real & epot,
                           real & coll_time, bool mw_flag, bool df_flag, const real vcirc,
                           const real ds)
{
    for (int i = 0; i < n ; i++)
        for (int k = 0; k < NDIM ; k++)
            acc[i][k] = jerk[i][k] = 0;
    epot = 0;
    const real VERY_LARGE_NUMBER = 1e300;
    real coll_time_q = VERY_LARGE_NUMBER;      // collision time to 4th power
    real coll_est_q;                           // collision time scale estimate
                                               // to 4th power (quartic)

    // G in units of kpc*km*km/s/s/Msun
    //const real G = 4.3e-6
    const real G = 4.4967e-6;              // G in kpc^3/Gyr^2/Msun
    const real pi = 3.141592653589793238462643383279;  // from www.math.com

    // Looping through the particles
    for (int i = 0; i < n ; i++){

  //  If the mass is non-zero (not a test particle)
  //  But must loop if we're using the Milky Way potential
  //  even if it's a test particles, b/w otherwise it won't "feel" the MW potential
	if (mass[i] != 0.0 || mw_flag == 1) {
	  // cerr << mass[i] << " " << i << endl;

        if (mw_flag) {

        //  MILKY WAY POTENTIAL

           // Parameters for the potential
           const real a = 6.5;
           const real b = 0.26;
	   const real rs = 0.7;    // c in paper
        //   const real c = 13.0;     // now ds, d in paper (12.0 kpc in paper), an input param
           const real q = 0.9;
           const real p = 1.0;

	   const real vcirc2 = vcirc*vcirc;          // v_sun^2
       //   const real vcirc2 = 225.0*225.0;          // v_sun^2
	   const real vh2 = 121.0*121.0;             // v_halo^2   114 km/s in paper
	   const real GM = 8.887*vcirc2;       // G*M_disk
	   const real GMs = 3.0*vcirc2;        // G*M_sphere
        
        //  ACCELERATION
        
	  // HALO
	  real x2 = pos[i][0]*pos[i][0];
	  real y2 = pos[i][1]*pos[i][1];
	  real z2 = pos[i][2]*pos[i][2];
	  real thx = 2.*vh2/(x2+y2/(p*p)+z2/(q*q)+ds*ds);
	//  real thx = 2.*vh2/(x2+y2/(p*p)+z2/(q*q)+c*c);
	  real thy = thx/(p*p);
	  real thz = thx/(q*q);
	  real phih = vh2*log(x2+y2/(p*p)+z2/(q*q)+ds*ds);
	//  real phih = vh2*log(x2+y2/(p*p)+z2/(q*q)+c*c);
	  real r2 = x2 + y2 + z2;
	  real thrad = 2.*vh2/(r2+ds*ds);
	//  real thrad = 2.*vh2/(r2+c*c);

          // SPHEROID
          real r = sqrt(r2);
          real rrs2 = (r+rs)*(r+rs);
          real tsrad = (GMs/rrs2)/r;
          real phis = -GMs/(r+rs);

          // DISK
          real rho = sqrt(x2 + y2);
          real Z = pos[i][2];
          real sqz2b2 = sqrt(Z*Z + b*b);
          real azb2 = (a + sqz2b2)*(a + sqz2b2);
          real tdr = GM/pow((rho*rho + azb2),1.5);
          real tdz = tdr*(a/sqz2b2 + 1.);
          real phim = -GM/sqrt(rho*rho+azb2);

          // Acceleration
          acc[i][0] += -(thx+tsrad+tdr)*pos[i][0];   // a_x
          acc[i][1] += -(thy+tsrad+tdr)*pos[i][1];   // a_y
          acc[i][2] += -(thz+tsrad+tdz)*pos[i][2];   // a_z

        // ADDING THE POTENTIAL
	epot+=mass[i]*(phih+phim+phis);


        // JERK

	  // HALO
 	  real term1 = pos[i][0]*vel[i][0] + pos[i][1]*vel[i][1]/(p*p) + pos[i][2]*vel[i][2]/(q*q);
          real thh = 4.*vh2*term1/pow(x2+y2/(p*p)+z2/(q*q)+ds*ds,2);
	//  real thh = 4.*vh2*term1/pow(x2+y2/(p*p)+z2/(q*q)+c*c,2);
 	  jerk[i][0] += -thx*vel[i][0] + thh*pos[i][0];
	  jerk[i][1] += -thx*vel[i][1]/(p*p) + thh*pos[i][1]/(p*p);
	  jerk[i][2] += -thx*vel[i][2]/(q*q) + thh*pos[i][2]/(q*q);

          // SPHEROID
          real rdotv = pos[i][0]*vel[i][0]+pos[i][1]*vel[i][1]+pos[i][2]*vel[i][2];
          real rrs3 = rrs2*(r+rs);
          real r3 = r*r*r;
          jerk[i][0] += -GMs*(vel[i][0]/r - rdotv*pos[i][0]/r3)/rrs2;
          jerk[i][0] += 2.*GMs*(rdotv*pos[i][0]/r2)/rrs3;
          jerk[i][1] += -GMs*(vel[i][1]/r - rdotv*pos[i][1]/r3)/rrs2;
          jerk[i][1] += 2.*GMs*(rdotv*pos[i][1]/r2)/rrs3;
          jerk[i][2] += -GMs*(vel[i][2]/r - rdotv*pos[i][2]/r3)/rrs2;
          jerk[i][2] += 2.*GMs*(rdotv*pos[i][2]/r2)/rrs3;

          // DISK
          real dend = sqrt(rho*rho + azb2);       // denominator
    	  real asq1 = (a/sqz2b2+1.);
          real vterm = pos[i][0]*vel[i][0] + pos[i][1]*vel[i][1] + asq1*pos[i][2]*vel[i][2];
          real zterm = (a*b*b)/(sqz2b2*sqz2b2*sqz2b2) + 1.;
          real dend3 = dend*dend*dend;
	  real dend5 = dend3*dend*dend;
          jerk[i][0] += -GM*(vel[i][0]/dend3 - 3.*pos[i][0]*vterm/dend5);
          jerk[i][1] += -GM*(vel[i][1]/dend3 - 3.*pos[i][1]*vterm/dend5);
          jerk[i][2] += -GM*(vel[i][2]*zterm/dend3 - 3.*pos[i][2]*asq1*vterm/dend5);
//          cerr << "zterm" << zterm << endl;
          //cerr << "jerkx" << jerk[i][0] << endl;
          //cerr << "jerky" << jerk[i][1] << endl;
          //cerr << "jerkz" << jerk[i][2] << endl;
         
	// CALCULATING COLLISION TIME

          // first collision time estimate, based on unaccelerated linear motion:

   	  real v2 = vel[i][0]*vel[i][0]+vel[i][1]*vel[i][1]+vel[i][2]*vel[i][2];
          coll_est_q = (r2*r2) / (v2*v2);
          if (coll_time_q > coll_est_q)
              coll_time_q = coll_est_q;

          // second collision time estimate, based on free fall:

          real da2 = 0;                                  // da2 becomes the 
          for (int k = 0; k < NDIM ; k++)                // square of the 
              da2 += acc[i][k] * acc[i][k];                      // pair-wise accel-

          coll_est_q = r2/da2;
          if (coll_time_q > coll_est_q)
              coll_time_q = coll_est_q;

          // cerr << v2 << " " << r2 << " " << da2 << " " << coll_time_q << endl;


          // DYNAMICAL FRICTION (mostly copied from orbit.f)
          if (df_flag) {
      
            // G in units of kpc*km*km/s/s/Msun
            //G = 4.3e-6
            //pi = 4.0*ATAN(1.0)
      
            real mfric = mass[i];
      
            // Calculate mass enclosed within r (=(1/r^2 G)(d phi /dr))
            real mr = 2.0*vh2*r3/(r2+ds*ds)/G;
      
            // Calculate speed - note: strictly should be updated so speed coincides with position
            real v = sqrt(v2);
      
            // Coulomb logarithm via James' fitting function
            real ot = 1.0/3.0;
            real rtide = r*pow((ot*mfric/mr),ot);
            real xlim = rtide;
            real coulog = log(r/rtide)+jfit(xlim)+0.135;
      
            // coulog is the Coulomb Logarithm
            // determined from Hashimoto 2003 
            // note: coulog becomes 0 if satellite radius is less than 1.4 epsilon
            // this is so we don't get dynamical acceleration at small radii
            // epsilon is the plummer scale length of the satellite
            real epsilon = soft[i];
            real coulog_hashi;
            if (r >= (1.4*epsilon)) {
              coulog_hashi = log(r / (1.4 * epsilon));
            } else {
              coulog_hashi = 0.;
            }      

            real df;
      
            if (coulog > 0.0) {
      
              // Local density (=(1/4 pi r^2)(d mr/dr))
              //real rho = (vh2/2.0d0/pi/G)*(r*r+3.0d0*c*c)/(r*r+c*c)/(r*r+c*c)
              real rho = (vh2/2.0/pi/G)*(r2+3.0*ds*ds)/(r2+ds*ds)/(r2+ds*ds);
      
              real sigma2 = vh2;
              real xx = v/sqrt(2.0*sigma2);  // capital X in B+T
      
              // B+T, eqn 7-18, df = dv/dt = acc without the vector part
              df = 4.0*pi*coulog*(G*G)*mfric*rho*(erf(xx)-2.0*xx*exp(-(xx*xx))/sqrt(pi))/(v*v*v);
       
            } else { 
              df = 0.0;
            }
       
            //dfx = -df*x(4)
            //dfy = -df*x(5)
            //dfz = -df*x(6)
      
            // Calculating the Acceleration
            acc[i][0] += -df*vel[i][0];
            acc[i][1] += -df*vel[i][1];
            acc[i][2] += -df*vel[i][2];
      
            // Calculating the Jerk
            if (coulog > 0.0) {
      
              real sigma2 = vh2;
              real xx = v/sqrt(2.0*sigma2);  // capital X in B+T
      
              real v3 = v2*v;
              real v4 = v2*v2;
              real vdotadv = ( vel[i][0]*acc[i][0]+vel[i][1]*acc[i][1]+vel[i][2]*acc[i][2] )/v;
              real j1 = -4.0*pi*coulog*(G*G)*mfric*rho;
              real j2 = 2.0*exp(-xx*xx)*vdotadv/sqrt(2.0*pi*sigma2);
              //real j2 = ( exp(-xx*xx)-exp(-2.0*xx)*(1.0+2.0*xx) )*2.0*vdotadv/sqrt(2.0*pi*sigma2);
              real j3 = erf(xx)-2.0*xx*exp(-(xx*xx))/sqrt(pi);
              jerk[i][0] += j1*(j2*vel[i][0]/v3 + j3*(acc[i][0]/v3-3.0*vel[i][0]*vdotadv/v4));
              jerk[i][1] += j1*(j2*vel[i][1]/v3 + j3*(acc[i][1]/v3-3.0*vel[i][1]*vdotadv/v4));
              jerk[i][2] += j1*(j2*vel[i][2]/v3 + j3*(acc[i][2]/v3-3.0*vel[i][2]*vdotadv/v4));
            }
      
          }    // end dynamical friction



        }    // end MW potential


  //  If the mass is non-zero (not a test particle)
	if (mass[i] != 0.0) {

  //  Looping through all pairs
        for (int j = i+1; j < n ; j++){            // rji[] is the vector from
    	    ///const real G = 4.4967e-6;              // G in kpc^3/Gyr^2/Msun, defined above

            real rji[NDIM];                        // particle i to particle j
            real vji[NDIM];                        // vji[] = d rji[] / d t
            for (int k = 0; k < NDIM ; k++){
                rji[k] = pos[j][k] - pos[i][k];
                vji[k] = vel[j][k] - vel[i][k];
            }
            real r2 = 0;                           // | rji |^2
            real v2 = 0;                           // | vji |^2
            real rv_r2 = 0;                        // ( rij . vij ) / | rji |^2
            real rv_re2 = 0;                       // ( rij . vij ) / ( |rji| * (rji + e) )
            for (int k = 0; k < NDIM ; k++){
                r2 += rji[k] * rji[k];
                v2 += vji[k] * vji[k];
                rv_r2 += rji[k] * vji[k];
            }
            rv_r2 /= r2;
            real r = sqrt(r2);                     // | rji |
            // rv_re2 /= r*(r+soft[i]);
            rv_re2 /= r*(r2+soft[i]*soft[i]);
            real r3 = r * r2;                      // | rji |^3
            // real re = r + soft[i];                 // |rji| + e,  e softening param
	    // real re3 = re*re*re;                   // ( |rji| + e )^3
            real sqr2e2 = sqrt(r2+soft[i]*soft[i]);
            real sqr2e2_3 = sqr2e2 * sqr2e2 * sqr2e2;

// add the {i,j} contribution to the total potential energy for the system:

            // epot -= G * mass[i] * mass[j] / r;
            // epot -= G * mass[i] * mass[j] / (r+soft[i]);
            epot -= G * mass[i] * mass[j] / sqr2e2;

// add the {j (i)} contribution to the {i (j)} values of acceleration and jerk:

            real da[3];                            // main terms in pairwise
            real dj[3];                            // acceleration and jerk
            for (int k = 0; k < NDIM ; k++){
	      // da[k] = rji[k] / r3;                           // see equations
	      // dj[k] = (vji[k] - 3 * rv_r2 * rji[k]) / r3;    // in the header
              // da[k] = rji[k] / re3;                           // see equations
              //  dj[k] = (vji[k] - 3 * rv_re2 * rji[k]) / re3;    // in the header
                da[k] = rji[k] / sqr2e2_3;       // a = -GM/(r^2 + e^2)^(3/2)
                dj[k] = (vji[k] - 3 * rv_re2 * rji[k]) / sqr2e2_3;    // in the header
            }

            for (int k = 0; k < NDIM ; k++){
                acc[i][k] += G * mass[j] * da[k];                 // using symmetry
                acc[j][k] -= G * mass[i] * da[k];                 // find pairwise
                jerk[i][k] += G * mass[j] * dj[k];                // acceleration
                jerk[j][k] -= G * mass[i] * dj[k];                // and jerk
            }

            //  If the mass is non-zero (not a test particle)
	    if (mass[j] != 0.0) {

              // first collision time estimate, based on unaccelerated linear motion:

              coll_est_q = (r2*r2) / (v2*v2);
              if (coll_time_q > coll_est_q)
                  coll_time_q = coll_est_q;

              // second collision time estimate, based on free fall:

              real da2 = 0;                                  // da2 becomes the 
              for (int k = 0; k < NDIM ; k++)                // square of the 
                  da2 += da[k] * da[k];                      // pair-wise accel-
              double mij = mass[i] + mass[j];                // eration between
              da2 *= G * G * mij * mij;                      // particles i and j

              coll_est_q = r2/da2;
              if (coll_time_q > coll_est_q)
                  coll_time_q = coll_est_q;

            }  // mass[j] != 0.0

        }  // for j

        }  // if mass[i] != 0

        }  // if mass[i] != 0 or MW on                           

    }     // for i                                          // from q for quartic back
    coll_time = sqrt(sqrt(coll_time_q));            // to linear collision time

    //cerr << coll_time << endl;
    //cerr << acc[0][0] << ' '<< acc[0][1] << ' ' << acc[0][2] << endl;
    //cerr << acc[1][0] << ' '<< acc[1][1] << ' ' << acc[1][2] << endl;
    //cerr << endl;
  
}                                             


//***********************************************************************
real jfit(real x)
{

real jfitv, num, den, fx;
real a,b,c,d,e,f,g,h,i,j;

// Fitting parameters.

     // PARAMETER (a=0.2021d0,b=0.04474d0,c=3.1621d0,
     //& d=0.5488d0,e=2.9965d0,f=0.7375d0,g=0.9217d0,
     //& h=0.2749d0,i=1.2583d0,j=2.3069d0)

      a=0.2021;
      b=0.04474;
      c=3.1621;
      d=0.5488;
      e=2.9965;
      f=0.7375;
      g=0.9217;
      h=0.2749;
      i=1.2583;
      j=2.3069;

      num = log(1.0+a*x)*(b*( pow(x,c) )+d*( pow(x,e) ));
      den = pow( 1.0 + f*( pow(x,g) ) + h*( pow(x,i) ) ,j);
      fx = log(1.0+x)-x/(1.0+x);

      jfitv = num/den;
      jfitv = jfitv/fx/fx;

      return jfitv;
}


/*-----------------------------------------------------------------------------
 *                                                                    \\   o
 *  end of file:  nbody_sh1.C                                         /\\'  O
 *                                                                   /\     |
 *=============================================================================
 */
