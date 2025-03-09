interface CardProps {
  title: string;
  content: string;
}

export function Card({ title, content }: CardProps) {
  return (
    <div className="overflow-hidden rounded-xl bg-white shadow-lg transition-all hover:shadow-xl">
      <div className="border-b border-gray-100 bg-gradient-to-r from-gray-50 to-white px-6 py-4">
        <h3 className="text-xl font-semibold text-gray-900">{title}</h3>
      </div>
      <div className="px-6 py-4">
        <p className="whitespace-pre-wrap text-gray-700 leading-relaxed">{content}</p>
      </div>
    </div>
  );
}