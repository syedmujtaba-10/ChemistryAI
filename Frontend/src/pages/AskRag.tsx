import { useState } from 'react';
import axios from 'axios';
import { Send } from 'lucide-react';
import { cn } from '../lib/utils';

interface Message {
  type: 'user' | 'assistant';
  content: string;
}

interface RagResponse {
  question: string;
  detailed_explanation: string;
}

export function AskRag() {
  const [messages, setMessages] = useState<Message[]>([]);
  const [input, setInput] = useState('');
  const [loading, setLoading] = useState(false);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    if (!input.trim() || loading) return;

    const question = input.trim();
    setInput('');
    setLoading(true);
    setMessages((prev) => [...prev, { type: 'user', content: question }]);

    try {
      const response = await axios.get<RagResponse>(
        `http://127.0.0.1:8000/chemistry/ask-rag?question=${encodeURIComponent(question)}`
      );
      setMessages((prev) => [
        ...prev,
        { type: 'assistant', content: response.data.detailed_explanation },
      ]);
    } catch (err) {
      setMessages((prev) => [
        ...prev,
        { type: 'assistant', content: 'Sorry, I encountered an error. Please try again.' },
      ]);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="flex h-full flex-col">
      <div className="flex-1 space-y-4 overflow-y-auto p-4">
        {messages.map((message, i) => (
          <div
            key={i}
            className={cn(
              'flex w-full',
              message.type === 'user' ? 'justify-end' : 'justify-start'
            )}
          >
            <div
              className={cn(
                'max-w-[80%] rounded-lg px-4 py-2',
                message.type === 'user'
                  ? 'bg-blue-600 text-white'
                  : 'bg-gray-100 text-gray-900'
              )}
            >
              <p className="whitespace-pre-wrap">{message.content}</p>
            </div>
          </div>
        ))}
        {loading && (
          <div className="flex justify-start">
            <div className="max-w-[80%] animate-pulse rounded-lg bg-gray-100 px-4 py-2">
              <div className="h-4 w-12 rounded bg-gray-300" />
            </div>
          </div>
        )}
      </div>

      <form onSubmit={handleSubmit} className="border-t p-4">
        <div className="flex gap-2">
          <input
            type="text"
            value={input}
            onChange={(e) => setInput(e.target.value)}
            placeholder="Ask a chemistry question..."
            className="flex-1 rounded-md border border-gray-300 px-3 py-2 shadow-sm focus:border-blue-500 focus:outline-none focus:ring-1 focus:ring-blue-500"
          />
          <button
            type="submit"
            disabled={loading}
            className="rounded-md bg-blue-600 p-2 text-white hover:bg-blue-700 focus:outline-none focus:ring-2 focus:ring-blue-500 focus:ring-offset-2 disabled:opacity-50"
          >
            <Send className="h-5 w-5" />
          </button>
        </div>
      </form>
    </div>
  );
}