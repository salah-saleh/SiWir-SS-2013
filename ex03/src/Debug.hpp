#ifndef DEBUG_HH
#define DEBUG_HH


//http://patrick.wagstrom.net/old/misc/dotfiles/dotbash_profile.html
#define BLACK         "\033[0;30m"
#define RED           "\033[0;31m"
#define GREEN         "\033[0;32m"
#define BROWN         "\033[0;33m"
#define BLUE          "\033[0;34m"
#define MAGENTA       "\033[0;35m"
#define CYAN          "\033[0;36m"
#define WHITE         "\033[0;37m"
#define LIGHTBLACK    "\033[1;30m"
#define LIGHTRED      "\033[1;31m"
#define LIGHTGREEN    "\033[1;32m"
#define YELLOW        "\033[1;33m"
#define LIGHTBLUE     "\033[1;34m"
#define LIGHTMAGENTA  "\033[1;35m"
#define LIGHTCYAN     "\033[1;36m"
#define LIGHTWHITE    "\033[1;37m"



#define ERROR(message)  ...
  
#ifdef DEBUG

// Asserts expr expression, if false exits with error
#define ASSERT(expr, message)  ...

#define INFO(message)  ...

#define WARN(message)  ...

#else
#define ASSERT(expr, message)
#define INFO(message)
#define WARN(message)                

#endif //DEBUG
#endif // DEBUG_HH
