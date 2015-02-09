/*
 * vtkdisplay.t
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

#include "globals.hpp"

#include <vtkVersion.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkImageImport.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkImageActor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkXMLImageDataWriter.h>

#include <boost/make_shared.hpp>

#include <cmath>

template<typename Type, typename PointerType>
void VTKGridDisplay<Type, PointerType>::write()
{
  imageImport->SetImportVoidPointer(this->field->getRawData());
  imageImport->Modified();
//  imageImport->GetInput()->Modified();
  imageImport->SetUpdateExtentToWholeExtent();
  imageImport->Update();
}

template<typename Type, typename PointerType>
void VTKGridDisplay<Type, PointerType>::close()
{

}

template<typename Type, typename PointerType>
void VTKGridDisplay<Type, PointerType>::init()
{
  schnek::SimpleDiagnostic<Type, PointerType>::init();
  int width  = this->field->getDims()[0];
  int height  = this->field->getDims()[1];

//  for (int i = this->field->getLo(0); i<=this->field->getHi(0); ++i)
//    for (int j = this->field->getLo(1); j<=this->field->getHi(1); ++j)
//    {
//      double x = i/(double)width;
//      double y = j/(double)height;
//      double r = sqrt(x*x + y*y);
//      (*(this->field))(i,j) = sin(12*r);
//    }

  imageImport = vtkSmartPointer<vtkImageImport>::New();

  imageImport->SetDataSpacing(1, 1, 1);
  imageImport->SetDataOrigin(0, 0, 0);
  imageImport->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
  imageImport->SetDataExtentToWholeExtent();
  imageImport->SetDataScalarType(TypeToVtk<typename Type::value_type>::value);
  imageImport->SetNumberOfScalarComponents(1);
  imageImport->SetImportVoidPointer(this->field->getRawData());
  imageImport->Update();

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

  Globals::instance().addInteractor(renderWindowInteractor);
  //renderWindowInteractor->Start();
}

