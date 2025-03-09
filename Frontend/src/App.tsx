import { BrowserRouter, Routes, Route, Outlet } from 'react-router-dom';
import { Login } from './pages/Login';
import { Register } from './pages/Register';
import { PredictProduct } from './pages/PredictProduct';
import { AskRag } from './pages/AskRag';
import { Flashcards } from './pages/Flashcards';
import { Navbar } from './components/Navbar';

// Layout that includes the Navbar and a main content area.
function MainLayout() {
  return (
    <div className="flex min-h-screen bg-gradient-to-br from-gray-50 via-white to-blue-50">
      <Navbar />
      <main className="ml-64 flex-1 p-8">
        <Outlet />
      </main>
    </div>
  );
}

function App() {
  return (
    <BrowserRouter>
      <Routes>
        {/* Routes without navbar */}
        <Route path="/login" element={<Login />} />
        <Route path="/register" element={<Register />} />

        {/* Routes with navbar */}
        <Route element={<MainLayout />}>
          <Route path="/predict-product" element={<PredictProduct />} />
          <Route path="/ask-rag" element={<AskRag />} />
          <Route path="/flashcards" element={<Flashcards />} />
        </Route>

        {/* Landing page is Login */}
        <Route path="/" element={<Login />} />
      </Routes>
    </BrowserRouter>
  );
}

export default App;
