#define CHECK( EXP )                     \
   do {                                                     \
      if( ! EXP ) {                             \
        Error( m_name.c_str(), "Failed to execute: %s at line %d", #EXP, __LINE__ );  \
        return EL::StatusCode::FAILURE;                    \
      }                                                     \
   } while( false )


