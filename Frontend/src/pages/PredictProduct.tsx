import { useState } from "react";
import axios from "axios";
import { supabase } from "../supabase"; // Adjust the path as needed
import { Loader2, FlaskConical } from "lucide-react";

interface PredictionResponse {
  reactants: string[];
  reactants_smiles: string; // e.g., "CCO.CC(=O)O"
  predicted_product_smiles: string; // e.g., "CC(=O)O.CCO>>CCOC(C)=O"
  human_readable_reaction: string;
  detailed_explanation: string;
  reaction_image: string; // base64 PNG data URL from backend
}

export function PredictProduct() {
  const [reactant1, setReactant1] = useState("");
  const [reactant2, setReactant2] = useState("");
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState<PredictionResponse | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [saveMessage, setSaveMessage] = useState<string>("");

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setLoading(true);
    setError(null);
    setResult(null);
    setSaveMessage("");

    const reactantsCombined = `${reactant1.trim()},${reactant2.trim()}`;

    try {
      const response = await axios.get<PredictionResponse>(
        `http://127.0.0.1:8000/chemistry/full-explanation?reactants=${encodeURIComponent(
          reactantsCombined
        )}`
      );
      setResult(response.data);
    } catch (err) {
      setError("Failed to predict reaction. Please try again.");
      console.error(err);
    } finally {
      setLoading(false);
    }
  };

  const handleSaveFlashcard = async () => {
    // Get the current user
    const {
      data: { user },
      error: userError,
    } = await supabase.auth.getUser();
    if (userError || !user) {
      setSaveMessage("You must be logged in to save a flashcard.");
      return;
    }

    // Prepare flashcard data; you can store the entire prediction result.
    const flashcardData = {
      user_id: user.id,
      data: result, // stores the entire prediction JSON
    };

    const { error: insertError } = await supabase
      .from("flashcards")
      .insert([flashcardData]);

    if (insertError) {
      console.error("Error saving flashcard:", insertError);
      setSaveMessage("Error saving flashcard.");
    } else {
      setSaveMessage("Flashcard saved successfully!");
    }
  };

  return (
    <div className="space-y-8">
      <div className="relative rounded-xl bg-gradient-to-br from-blue-50 to-indigo-50 p-8 shadow-lg">
        <div className="absolute right-4 top-4">
          <FlaskConical className="h-12 w-12 text-indigo-300" />
        </div>
        <h1 className="mb-6 text-3xl font-bold text-gray-900">
          Chemical Reaction Predictor
        </h1>
        <form onSubmit={handleSubmit} className="space-y-6">
          <div className="grid gap-6 md:grid-cols-2">
            <div className="relative">
              <label
                htmlFor="reactant1"
                className="mb-2 block text-sm font-medium text-gray-700"
              >
                First Reactant
              </label>
              <input
                type="text"
                id="reactant1"
                value={reactant1}
                onChange={(e) => setReactant1(e.target.value)}
                placeholder="e.g., ethanol"
                required
                className="block w-full rounded-lg border border-gray-300 bg-white p-2.5"
              />
            </div>
            <div className="relative">
              <label
                htmlFor="reactant2"
                className="mb-2 block text-sm font-medium text-gray-700"
              >
                Second Reactant
              </label>
              <input
                type="text"
                id="reactant2"
                value={reactant2}
                onChange={(e) => setReactant2(e.target.value)}
                placeholder="e.g., acetic acid"
                required
                className="block w-full rounded-lg border border-gray-300 bg-white p-2.5"
              />
            </div>
          </div>
          <button
            type="submit"
            disabled={loading}
            className="w-full rounded-lg bg-indigo-500 p-2 text-white"
          >
            {loading ? (
              <>
                <Loader2 className="inline h-5 w-5 animate-spin" />
                <span> Predicting...</span>
              </>
            ) : (
              <>
                <span>Predict Reaction</span>
              </>
            )}
          </button>
        </form>
      </div>

      {error && (
        <div className="rounded-lg bg-red-50 p-4 text-red-700 shadow-md">
          <p className="font-medium">{error}</p>
        </div>
      )}

      {result && (
        <div className="space-y-6">
          <div className="rounded-lg bg-white p-4 shadow-md">
            <h2 className="text-xl font-semibold">Prediction Details</h2>
            <p>
              <strong>Reactants:</strong> {result.reactants.join(" + ")}
            </p>
            <p>
              <strong>Reactants (SMILES):</strong> {result.reactants_smiles}
            </p>
            <p>
              <strong>Predicted Product (SMILES):</strong>{" "}
              {result.predicted_product_smiles}
            </p>
          </div>

          <div className="rounded-lg bg-white p-4 shadow-md">
            <h2 className="text-xl font-semibold">Human-Readable Reaction</h2>
            <p>{result.human_readable_reaction}</p>
          </div>

          <div className="rounded-lg bg-white p-4 shadow-md">
            <h2 className="text-xl font-semibold">Detailed Explanation</h2>
            <div className="whitespace-pre-wrap prose">
              {result.detailed_explanation}
            </div>
          </div>

          <div className="rounded-lg bg-gray-100 p-4 shadow-md">
            <h2 className="mb-3 text-xl font-semibold text-gray-900">
              Reaction Diagram
            </h2>
            {result.reaction_image ? (
              <img
                src={result.reaction_image}
                alt="Reaction Diagram"
                className="w-full h-auto border rounded-lg"
              />
            ) : (
              <p>Unable to generate reaction image.</p>
            )}
          </div>

          <div className="flex justify-end">
            <button
              onClick={handleSaveFlashcard}
              className="rounded-lg bg-green-500 px-4 py-2 text-white hover:bg-green-600"
            >
              Save Flashcard
            </button>
            {saveMessage && (
              <p className="ml-4 self-center text-sm text-gray-700">
                {saveMessage}
              </p>
            )}
          </div>
        </div>
      )}
    </div>
  );
}
