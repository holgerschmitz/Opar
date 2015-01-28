/*
 * vtkdisplay.cpp
 *
 * Created on: 20 Jan 2015
 * Author: Holger Schmitz
 * Email: holger@notjustphysics.com
 *
 * Copyright 2012 Holger Schmitz
 *
 * This file is part of OPar.
 *
 * OPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * OPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with OPar.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "vtkdisplay.hpp"

#include <vtkVersion.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkImageImport.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkImageActor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkXMLImageDataWriter.h>

void VTKGridDisplay::write()
{

}

void VTKGridDisplay::close()
{

}

void VTKGridDisplay::init()
{

}


int main(int, char *[])
{
  // Create a c-style image
  const int width = 4;
  const int height = 4;

  unsigned char cImage[width*height];
  unsigned char value = 0;
  for(unsigned int row = 0; row < height; ++row)
    {
    for(unsigned int col = 0; col < width; ++col)
      {
      cImage[row * width + col] = value;
      value += 10;
      }
    }

  // Convert the c-style image to a vtkImageData
  vtkSmartPointer<vtkImageImport> imageImport =
    vtkSmartPointer<vtkImageImport>::New();
  imageImport->SetDataSpacing(1, 1, 1);
  imageImport->SetDataOrigin(0, 0, 0);
  imageImport->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
  imageImport->SetDataExtentToWholeExtent();
  imageImport->SetDataScalarTypeToUnsignedChar();
  imageImport->SetNumberOfScalarComponents(1);
  imageImport->SetImportVoidPointer(cImage);
  imageImport->Update();

  // Create an actor
  vtkSmartPointer<vtkImageActor> actor =
    vtkSmartPointer<vtkImageActor>::New();
#if VTK_MAJOR_VERSION <= 5
  actor->SetInput(imageImport->GetOutput());
#else
  actor->SetInputData(imageImport->GetOutput());
#endif

  // Setup renderer
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  renderer->AddActor(actor);
  renderer->ResetCamera();

  // Setup render window
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);

  // Setup render window interactor
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  vtkSmartPointer<vtkInteractorStyleImage> style =
    vtkSmartPointer<vtkInteractorStyleImage>::New();

  renderWindowInteractor->SetInteractorStyle(style);

  // Render and start interaction
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderWindowInteractor->Initialize();

  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
