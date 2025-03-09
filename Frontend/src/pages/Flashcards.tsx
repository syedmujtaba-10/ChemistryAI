import { useState, useEffect } from "react";
import { supabase } from "../supabase"; // Adjust path as needed
import { Loader2, Trash } from "lucide-react";

interface PredictionResponse {
  reactants: string[];
  reactants_smiles: string;
  predicted_product_smiles: string;
  human_readable_reaction: string;
  detailed_explanation: string;
  reaction_image: string;
}

interface Flashcard {
  id: string;
  user_id: string;
  data: PredictionResponse;
  created_at: string;
}

export function Flashcards() {
  const [flashcards, setFlashcards] = useState<Flashcard[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // Fetch flashcards for the current user
  const fetchFlashcards = async () => {
    setLoading(true);
    setError(null);
    try {
      const {
        data: { user },
        error: userError,
      } = await supabase.auth.getUser();
      if (userError || !user) {
        throw new Error("User not logged in");
      }
      const userId = user.id;
      const { data, error } = await supabase
        .from("flashcards")
        .select("*")
        .eq("user_id", userId)
        .order("created_at", { ascending: false });
      if (error) {
        throw error;
      }
      setFlashcards(data as Flashcard[]);
    } catch (err: any) {
      setError(err.message || "Error fetching flashcards");
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    fetchFlashcards();
  }, []);

  // Delete a flashcard by its id
  const handleDelete = async (flashcardId: string) => {
    try {
      const { error } = await supabase
        .from("flashcards")
        .delete()
        .eq("id", flashcardId);
      if (error) throw error;
      // Refresh flashcards after deletion
      fetchFlashcards();
    } catch (err: any) {
      console.error("Error deleting flashcard:", err);
    }
  };

  return (
    <div className="p-8">
      <h1 className="mb-4 text-3xl font-bold">Your Flashcards</h1>
      {loading ? (
        <div className="flex items-center">
          <Loader2 className="h-6 w-6 animate-spin" />
          <span className="ml-2">Loading flashcards...</span>
        </div>
      ) : error ? (
        <p className="text-red-500">{error}</p>
      ) : flashcards.length === 0 ? (
        <p>No flashcards found.</p>
      ) : (
        <div className="space-y-4">
          {flashcards.map((card) => (
            <div key={card.id} className="rounded-lg bg-white p-4 shadow-md">
              <div className="flex items-center justify-between">
                <h2 className="text-xl font-semibold">Flashcard</h2>
                <button
                  onClick={() => handleDelete(card.id)}
                  className="rounded bg-red-500 p-2 text-white hover:bg-red-600"
                >
                  <Trash className="h-5 w-5" />
                </button>
              </div>
              <div className="mt-2 space-y-1">
                <p>
                  <strong>Reactants:</strong> {card.data.reactants.join(" + ")}
                </p>
                <p>
                  <strong>Reactants (SMILES):</strong> {card.data.reactants_smiles}
                </p>
                <p>
                  <strong>Predicted Product (SMILES):</strong>{" "}
                  {card.data.predicted_product_smiles}
                </p>
                <p>
                  <strong>Human-Readable Reaction:</strong>{" "}
                  {card.data.human_readable_reaction}
                </p>
                <div>
                  <strong>Detailed Explanation:</strong>
                  <div className="whitespace-pre-wrap mt-1 border-t pt-1">
                    {card.data.detailed_explanation}
                  </div>
                </div>
                {card.data.reaction_image && (
                  <div className="mt-4">
                    <img
                      src={card.data.reaction_image}
                      alt="Reaction Diagram"
                      className="w-full h-auto border rounded-lg"
                    />
                  </div>
                )}
              </div>
            </div>
          ))}
        </div>
      )}
    </div>
  );
}
