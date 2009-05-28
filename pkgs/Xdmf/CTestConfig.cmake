## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## # The following are required to uses Dart and the Cdash dashboard
##   ENABLE_TESTING()
##   INCLUDE(CTest)
set(CTEST_PROJECT_NAME "Xdmf")
set(CTEST_NIGHTLY_START_TIME "21:00:00 EST")

set(CTEST_DROP_METHOD "https")
set(CTEST_DROP_SITE "hsai-web-origin.arl.army.mil")
set(CTEST_DROP_LOCATION "/hsai/CDash/submit.php?project=Xdmf")
set(CTEST_DROP_SITE_CDASH TRUE)
set(CTEST_CURL_OPTIONS "CURLOPT_SSL_VERIFYPEER_OFF;CURLOPT_SSL_VERIFYHOST_OFF")
