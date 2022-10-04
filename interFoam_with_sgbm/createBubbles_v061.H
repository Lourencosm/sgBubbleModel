//Info << "---->createBubbles_v061.H: autoPtr to bubble model: start" << endl;

// Construct bubble model
autoPtr<bubbleModel> bubbleModel_
(
    bubbleModel::New(mesh)
);

//Info << "---->createBubbles_v061.H: autoPtr to bubble model: done" << endl;
