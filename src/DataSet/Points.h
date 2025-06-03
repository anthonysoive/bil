#ifndef POINTS_H
#define POINTS_H


/* Forward declarations */
struct Points_t; //typedef struct Points_t       Points_t ;
struct DataFile_t;
struct Mesh_t;
struct Point_t;


extern Points_t*  (Points_New)     (const int) ;
extern Points_t*  (Points_Create)  (DataFile_t*,Mesh_t*) ;
extern void       (Points_Delete)  (void*) ;



#define Points_GetNbOfPoints(PTS)    ((PTS)->n_points)
#define Points_GetPoint(PTS)         ((PTS)->point)


struct Points_t {
  int n_points ;              /* nb of points */
  Point_t*  point ;           /* Point */
} ;


#endif
