#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include "imgui/imgui.h"
#include <iostream>
#include "mxuSurfaceSimplification.h"
#include <chrono>

namespace GLOBAL_VARS {
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	int targetNumFaces = 1000;
}

Eigen::MatrixXd V2;
Eigen::MatrixXi F2;

int main(int argc, char* argv[])
{


	// Load a mesh in OFF format
	igl::readOFF("../data/bunny.off", GLOBAL_VARS::V, GLOBAL_VARS::F);
	//igl::readOFF(TUTORIAL_SHARED_PATH "/lion.off", V2, F2);

	int numFaces = GLOBAL_VARS::F.rows();
	int numVertices = GLOBAL_VARS::V.rows();

	// Init the viewer
	igl::opengl::glfw::Viewer viewer;

	// Attach a menu plugin
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	// Customize the menu
	double doubleVariable = 0.1f; // Shared between two menus

	// Add content to the default menu window
	menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu.draw_viewer_menu();

		// Add new group
		if (ImGui::CollapsingHeader("Surface Simplification", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::Text(std::string(std::string("number of faces: ") + std::to_string(numFaces)).c_str());
			ImGui::Text(std::string(std::string("number of vertices: ") + std::to_string(numVertices)).c_str());
			ImGui::InputInt("targetNumFaces", &GLOBAL_VARS::targetNumFaces);
			if (ImGui::Button("Reset Mesh")) {
				igl::readOFF("../data/bunny.off", GLOBAL_VARS::V, GLOBAL_VARS::F);
				numVertices = GLOBAL_VARS::V.rows();
				numFaces = GLOBAL_VARS::F.rows();

				viewer.data().clear();
				viewer.data().set_mesh(GLOBAL_VARS::V, GLOBAL_VARS::F);
			}
			if (ImGui::Button("Iterative Edge Contraction")) {
				IterativeEdgeContraction(GLOBAL_VARS::V, GLOBAL_VARS::F, GLOBAL_VARS::targetNumFaces);
				numVertices = GLOBAL_VARS::V.rows();
				numFaces = GLOBAL_VARS::F.rows();

				viewer.data().clear();
				viewer.data().set_mesh(GLOBAL_VARS::V, GLOBAL_VARS::F);
			}
			if (ImGui::Button("mxu Iterative Edge Contraction")) {
				MxuIterativeEdgeContraction(GLOBAL_VARS::V, GLOBAL_VARS::F, GLOBAL_VARS::targetNumFaces);
				numVertices = GLOBAL_VARS::V.rows();
				numFaces = GLOBAL_VARS::F.rows();

				viewer.data().clear();
				viewer.data().set_mesh(GLOBAL_VARS::V, GLOBAL_VARS::F);
			}
		}
	};

	// Draw additional windows
	menu.callback_draw_custom_window = [&]()
	{
		// Define next window position + size
		ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiSetCond_FirstUseEver);
		ImGui::Begin(
			"New Window", nullptr,
			ImGuiWindowFlags_NoSavedSettings
		);

		//// Expose the same variable directly ...
		//ImGui::PushItemWidth(-80);
		//ImGui::DragScalar("double", ImGuiDataType_Double, &doubleVariable, 0.1, 0, 0, "%.4f");
		//ImGui::PopItemWidth();

		static std::string str = "bunny";
		ImGui::InputText("Name", str);

		ImGui::End();
	};

	// Plot the mesh
	viewer.data().set_mesh(GLOBAL_VARS::V, GLOBAL_VARS::F);
	viewer.launch();
}
