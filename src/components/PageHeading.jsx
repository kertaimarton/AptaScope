export default function PageHeading({ children }) {
  return (
    <div>
      <h1 className="text-2xl font-semibold">{children}</h1>
      <div className="h-[2px] w-10 bg-accent mt-2" />
    </div>
  );
}
