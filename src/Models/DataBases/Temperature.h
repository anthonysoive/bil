#ifndef TEMPERATURE_H
#define TEMPERATURE_H


#define Temperature_0C                  (273.15)
#define Temperature_20C                 (293.15)
#define Temperature_25C                 (298.15)

#define Temperature_DefaultValue        (Temperature_20C)


#define Temperature_GetRoomValue(TEMP)     ((TEMP)->roomvalue)


#define Temperature_SetRoomTemperature(TEMP,T) \
        do {Temperature_GetRoomValue(TEMP) = (T) ;} while(0)


struct Temperature_t {
  double roomvalue ;          /* Room temperature */
  
  Temperature_t() {roomvalue = Temperature_DefaultValue;}
  ~Temperature_t(){};
};



inline Temperature_t* (Temperature_Create)(void){
  Temperature_t* temperature = new Temperature_t();
  
  return(temperature) ;
}



inline void (Temperature_Delete)(void* self){
  Temperature_t* temperature = (Temperature_t*) self;
  
  delete(temperature);
}

#endif
