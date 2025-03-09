import { useState, useEffect } from "react";
import { NavLink } from "react-router-dom";
import { supabase } from "../supabase"; // Adjust the path as needed
import { BeakerIcon, BrainCircuitIcon, Car as Cards } from "lucide-react";
import { cn } from "../lib/utils";

export function Navbar() {
  const links = [
    { to: "/predict-product", icon: BeakerIcon, label: "Predict Product" },
    { to: "/ask-rag", icon: BrainCircuitIcon, label: "Ask RAG" },
    { to: "/flashcards", icon: Cards, label: "Flashcards" },
  ];

  const [user, setUser] = useState<any>(null);

  useEffect(() => {
    // Get the current session on mount
    const getUser = async () => {
      const {
        data: { session },
      } = await supabase.auth.getSession();
      setUser(session?.user ?? null);
    };

    getUser();

    // Listen for auth state changes
    const {
      data: { subscription },
    } = supabase.auth.onAuthStateChange((event, session) => {
      setUser(session?.user ?? null);
    });

    return () => {
      subscription.unsubscribe();
    };
  }, []);

  return (
    <nav className="fixed left-0 top-0 flex h-full w-64 flex-col bg-gray-900 p-4">
      <div className="flex flex-col gap-2">
        {links.map(({ to, icon: Icon, label }) => (
          <NavLink
            key={to}
            to={to}
            className={({ isActive }) =>
              cn(
                "flex items-center gap-3 rounded-lg px-3 py-2 text-gray-300 transition-colors hover:bg-gray-800 hover:text-white",
                isActive && "bg-gray-800 text-white"
              )
            }
          >
            <Icon className="h-5 w-5" />
            <span>{label}</span>
          </NavLink>
        ))}
      </div>

      <div className="mt-auto border-t border-gray-700 pt-4">
        {user ? (
          <div className="flex flex-col gap-2">
            <span className="text-sm text-gray-300">Hello, {user.email}</span>
            <button
              className="flex items-center gap-2 rounded-lg bg-red-600 px-3 py-2 text-white hover:bg-red-700"
              onClick={async () => {
                await supabase.auth.signOut();
                setUser(null);
              }}
            >
              Logout
            </button>
          </div>
        ) : (
          <NavLink
            to="/login"
            className="flex items-center gap-2 rounded-lg px-3 py-2 text-gray-300 transition-colors hover:bg-gray-800 hover:text-white"
          >
            Login
          </NavLink>
        )}
      </div>
    </nav>
  );
}
