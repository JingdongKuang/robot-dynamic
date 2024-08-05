# CMakeList.txt: vsCollisionDetection 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.8)

# 如果支持，请为 MSVC 编译器启用热重载。
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("vsCollisionDetection")

# 将源代码添加到此项目的可执行文件。
add_executable (vsCollisionDetection "DH_dynamics.cpp" "DH_dynamics.h" "Chebyshef_40hz_Order3.cpp" "Chebyshef_40hz_Order3.h" "rtwtypes.h" "mian_dynamics.cpp" "mian_dynamics.h" "POE_kinematics.h" "POE_kinematics.cpp")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET vsCollisionDetection PROPERTY CXX_STANDARD 20)
endif()

# 指定Eigen库的路径
set(EIGEN3_INCLUDE_DIR "F:/externalLib/cpp/eigen-3.4.0")#此处不能用\,只能用/

# 将Eigen库链接到项目中
target_include_directories(vsCollisionDetection PRIVATE ${EIGEN3_INCLUDE_DIR})