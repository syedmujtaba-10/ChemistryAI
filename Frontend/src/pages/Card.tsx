import React from 'react';

interface CardProps {
  title: string;
  content: React.ReactNode; // Allows strings, JSX, and HTML elements
}

export function Card({ title, content }: CardProps) {
  return (
    <div className="bg-white rounded-lg shadow-lg p-6 border border-gray-200">
      <h2 className="text-xl font-semibold text-gray-800 mb-3">{title}</h2>
      <div className="text-gray-700">{content}</div>
    </div>
  );
}
