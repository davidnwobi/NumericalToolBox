if(EXISTS "C:/Users/user/OneDrive/Desktop/object_oriented_numerical_anaylsis/NumericalToolbox/build/tests/NumLibTests[1]_tests.cmake")
  include("C:/Users/user/OneDrive/Desktop/object_oriented_numerical_anaylsis/NumericalToolbox/build/tests/NumLibTests[1]_tests.cmake")
else()
  add_test(NumLibTests_NOT_BUILT NumLibTests_NOT_BUILT)
endif()
