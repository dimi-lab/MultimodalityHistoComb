function createForm() {
  // Create a new Google Form
  const form = FormApp.create("Virtual H&E Quiz");

  // Set form title and description
  form.setTitle("Image Quiz");
  form.setDescription("Each question includes an image and three multiple-choice options.");

  // Define choices for each question
  const choices = ['Real', 'Not sure', 'Fake'];

  // Google Drive folder ID containing images
  const folderId = '1w33i934Gvu-pbt53Gsd9bWHKl1JFejFw';  // Replace with your actual folder ID
  const folder = DriveApp.getFolderById(folderId);

  // Get all image files in the folder
  const files = folder.getFilesByType(MimeType.PNG);

  var cnt = 0;
  // Loop through the files and create a question for each image
  while (files.hasNext()) {
    const file = files.next();

    // Add the image as an ImageItem
    const img = file.getBlob(); // Retrieve the image directly from Drive
    const imageItem = form.addImageItem();
    imageItem.setImage(img).setTitle("Evaluate Image " + cnt);

    // Add a multiple-choice question after the image
    const questionItem = form.addMultipleChoiceItem();
    questionItem.setTitle("Is image " + cnt + " real or fake?")
                .setChoices(choices.map(choice => questionItem.createChoice(choice)));
    Logger.log(file.getName() + "\n");
    cnt += 1;
  }

  Logger.log("Form created: " + form.getEditUrl());
}
